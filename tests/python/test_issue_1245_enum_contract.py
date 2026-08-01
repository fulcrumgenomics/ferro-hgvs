"""Issue #1245: the exported enums must honour the contract their stubs declare.

All nine enums are pyo3 simple enums. pyo3's ``eq`` option generates ``__eq__``,
and a Python type that defines ``__eq__`` without ``__hash__`` is **unhashable**
— so eight of the nine crashed in idiomatic code like ``Counter(...)``,
``{x}`` or a dict keyed by the enum. ``Axis`` was the lone exception because it
alone carried ``frozen, hash``, which is what showed the omission was accidental
rather than designed.

The same block adds ordering (``ord``), ``.name`` / ``.value`` accessors, and
by-value construction, so the runtime now matches what the stubs promise.

Two further gaps came out of review, both about a member and the integer it
equals drifting apart:

* ``eq_int`` makes members compare equal to plain ints, but pyo3's ``hash``
  option hashes the *variant*, not its discriminant — so ``Axis.Genomic == 0``
  while ``hash(Axis.Genomic) != hash(0)``, and the two landed in different dict
  buckets. Each enum now defines ``__hash__`` as its discriminant.
* ``Cls(value)`` returned a fresh equal object rather than the interned member,
  so ``is`` — which the standard library documents as the way to compare enum
  members — silently disagreed with ``==``. The constructor now returns the
  interned member.

Two ``enum.Enum`` behaviours remain out of reach and the stubs no longer claim
them: iterating the class (``list(X)``, ``__members__``) needs a metaclass, and
``isinstance(x, enum.Enum)`` needs real ``enum.Enum`` inheritance. pyo3 provides
neither for a ``#[pyclass]`` enum.

These tests are written against every exported enum rather than a sample, since
the defect was precisely that one of the nine had been treated differently.
"""

import collections
import enum
import itertools
import pickle
from collections.abc import Callable

import pytest

import ferro_hgvs

ENUM_NAMES = [
    "Axis",
    "Consequence",
    "Impact",
    "ErrorMode",
    "ErrorType",
    "ErrorOverride",
    "EquivalenceLevel",
    "GenomeBuild",
    "Strand",
]

ENUMS = [getattr(ferro_hgvs, name) for name in ENUM_NAMES]


def members(cls: type) -> list[object]:
    """Every member bound as an attribute of ``cls``."""
    found = [
        getattr(cls, attribute)
        for attribute in dir(cls)
        if not attribute.startswith("__") and isinstance(getattr(cls, attribute, None), cls)
    ]
    assert found, f"{cls.__name__} exposes no members"
    return found


@pytest.fixture(params=ENUMS, ids=ENUM_NAMES)
def enum_class(request: pytest.FixtureRequest) -> type:
    return request.param


# ---------------------------------------------------------------------------
# The crashes.
# ---------------------------------------------------------------------------


def test_members_are_hashable(enum_class: type) -> None:
    for member in members(enum_class):
        hash(member)


def test_members_work_as_set_and_dict_keys(enum_class: type) -> None:
    all_members = members(enum_class)
    assert len(set(all_members)) == len(all_members)
    assert len({member: index for index, member in enumerate(all_members)}) == len(all_members)


def test_members_can_be_counted(enum_class: type) -> None:
    # The reported failure: Counter over a stream of enum values.
    first = members(enum_class)[0]
    assert collections.Counter([first, first])[first] == 2


def test_members_can_be_grouped(enum_class: type) -> None:
    first = members(enum_class)[0]
    grouped = [(key, len(list(group))) for key, group in itertools.groupby([first, first])]
    assert grouped == [(first, 2)]


def test_members_are_orderable(enum_class: type) -> None:
    all_members = members(enum_class)
    # Sorting must agree with the integer values, as IntEnum ordering would.
    assert [member.value for member in sorted(all_members)] == sorted(
        member.value for member in all_members
    )


def test_ordering_accepts_ints_and_rejects_foreign_operands(enum_class: type) -> None:
    # The operand set the stubs declare for `__lt__`/`__le__`/`__gt__`/`__ge__`
    # (`Self | int`). It was `Any`, which type-checked every one of the
    # rejections below. A self-bound TypeVar would go too far the other way and
    # reject the int comparison that `eq_int` + `ord` genuinely support.
    member = members(enum_class)[0]

    # An int operand works, and orders by the member's value.
    assert (member < member.value + 1) is True
    assert (member > member.value - 1) is True

    # A member of a *different* native enum does not — ordering is per-enum.
    other_enum = next(cls for cls in ENUMS if cls is not enum_class)
    with pytest.raises(TypeError):
        member < members(other_enum)[0]  # noqa: B015 - raising is the assertion

    # Nor does an unrelated type.
    with pytest.raises(TypeError):
        member < "not a member"  # noqa: B015 - raising is the assertion

    # Equality stays total across enums — `False`, not `TypeError` — which is
    # why `__eq__` keeps `object` rather than narrowing with the ordering ops.
    assert (member == members(other_enum)[0]) is False


def test_equal_members_hash_equally(enum_class: type) -> None:
    # The hash/eq contract: a == b implies hash(a) == hash(b).
    for member in members(enum_class):
        rebuilt = enum_class(member.value)
        assert rebuilt == member
        assert hash(rebuilt) == hash(member)


def test_members_hash_like_the_int_they_equal(enum_class: type) -> None:
    # `eq_int` extends equality across the type boundary — a member equals a
    # plain int — so the hash contract has to reach across it too. pyo3's
    # `hash` option does not: it delegates to the Rust `Hash` derive, which
    # hashes the variant rather than its discriminant. That left every member
    # equal to an int it hashed differently from.
    for member in members(enum_class):
        assert member == member.value
        assert hash(member) == hash(member.value), (
            f"{enum_class.__name__}.{member.name} == {member.value} but hashes differently"
        )


def test_members_and_their_ints_are_interchangeable_as_keys(enum_class: type) -> None:
    # The practical consequence, and how the defect actually surfaced: a
    # container keyed by the int did not answer to the member, and vice versa,
    # even though `==` said they were the same key.
    for member in members(enum_class):
        assert member in {member.value: "by int"}
        assert member.value in {member: "by member"}
        assert member.value in {member}
        assert member in {member.value}
        # A member and its int must collapse to one entry, never two.
        assert len({member, member.value}) == 1


# ---------------------------------------------------------------------------
# The accessors.
# ---------------------------------------------------------------------------


def test_members_expose_name(enum_class: type) -> None:
    for member in members(enum_class):
        assert isinstance(member.name, str)
        assert member.name
        # `.name` must be the attribute the member is actually bound to.
        assert getattr(enum_class, member.name) == member


def test_members_expose_value(enum_class: type) -> None:
    for member in members(enum_class):
        assert isinstance(member.value, int)
        assert member.value == int(member)


def test_values_are_distinct(enum_class: type) -> None:
    values = [member.value for member in members(enum_class)]
    assert len(set(values)) == len(values)


# ---------------------------------------------------------------------------
# By-value construction.
# ---------------------------------------------------------------------------


def test_every_member_round_trips_through_its_value(enum_class: type) -> None:
    # This is the guard on the generated constructors: a wrong or missing arm
    # for any of the 94 variants across the nine enums shows up right here.
    for member in members(enum_class):
        assert enum_class(member.value) == member
        assert enum_class(member.value).name == member.name


def test_construction_returns_the_interned_member(enum_class: type) -> None:
    # The standard library documents identity as *the* way to compare enum
    # members (`Color.RED is Color.RED`), so a constructor handing back a fresh
    # equal object is a footgun for anyone carrying that habit over: `==` would
    # work and `is` would silently not. Returning the interned member removes
    # the discrepancy rather than documenting around it.
    for member in members(enum_class):
        assert enum_class(member.value) is member
        # Repeated construction is the same object every time, not merely equal.
        assert enum_class(member.value) is enum_class(member.value)


def test_members_returned_from_the_library_are_also_interned(enum_class: type) -> None:
    # Interning has to hold for values crossing the Rust boundary too, not just
    # the constructor — that is where users actually get their members from.
    for member in members(enum_class):
        assert getattr(enum_class, member.name) is member


def test_construction_rejects_an_unknown_value(enum_class: type) -> None:
    unknown = max(member.value for member in members(enum_class)) + 1
    with pytest.raises(ValueError):
        enum_class(unknown)


def test_error_type_rejects_the_retired_discriminant() -> None:
    # `max(value) + 1` above only covers the value past the end, so it cannot
    # see a hole in the middle. Discriminant 5 was `PositionZero` (W4002),
    # retired in issue #269; the numeric value is skipped deliberately so the
    # surviving variants keep their stable repr-C discriminants. Nothing else
    # pins the gap, so a regression that reintroduced 5 would pass every other
    # test here.
    assert 5 not in {member.value for member in members(ferro_hgvs.ErrorType)}
    with pytest.raises(ValueError):
        ferro_hgvs.ErrorType(5)


# ---------------------------------------------------------------------------
# Behaviours that already worked and must keep working.
# ---------------------------------------------------------------------------


def test_members_compare_equal_to_their_int(enum_class: type) -> None:
    for member in members(enum_class):
        assert member == member.value


def test_members_render(enum_class: type) -> None:
    for member in members(enum_class):
        assert str(member)
        assert repr(member)


# ---------------------------------------------------------------------------
# The documented residue: what these are still not.
# ---------------------------------------------------------------------------


def test_members_are_not_enum_enum_instances(enum_class: type) -> None:
    # Pinned deliberately. pyo3 cannot make a #[pyclass] enum inherit from
    # enum.Enum, so the stubs no longer declare IntEnum. If a future pyo3 makes
    # this possible, this test is the reminder to revisit the stubs.
    assert not isinstance(members(enum_class)[0], enum.Enum)


def test_the_class_is_not_iterable(enum_class: type) -> None:
    # Likewise: iterating the class needs a metaclass pyo3 does not provide.
    with pytest.raises(TypeError):
        list(enum_class)


# ---------------------------------------------------------------------------
# Strand.Unknown was missing from the stub entirely (sibling defect).
# ---------------------------------------------------------------------------


def test_strand_exposes_all_three_members() -> None:
    assert {member.name for member in members(ferro_hgvs.Strand)} == {
        "Plus",
        "Minus",
        "Unknown",
    }


def test_stub_lists_exactly_the_runtime_members(
    enum_class: type, stub_enum_members: Callable[[str], dict[str, int]]
) -> None:
    """Every native enum's stub must match its runtime members.

    An equivalent guard already existed, but only for ``ErrorType`` (added for
    issue #630). That is precisely why ``Strand.Unknown`` could exist at runtime
    while the stub listed only ``Plus`` and ``Minus`` — nothing checked the
    other eight. Applying it to all of them closes the gap for good.
    """
    stub = stub_enum_members(enum_class.__name__)
    runtime = {member.name: member.value for member in members(enum_class)}
    assert stub == runtime


def test_members_are_not_picklable(enum_class: type) -> None:
    # Pinned, not skipped. Hashability often travels with pickling in consumer
    # code (caches keyed by enum, multiprocessing), so which way this behaves is
    # worth recording — but an `except Exception: pytest.skip` records nothing:
    # it reports success whether pickling works or not, so it can never fail.
    #
    # A pyo3 `#[pyclass]` enum has no `__reduce__`, so `pickle.dumps` raises
    # `TypeError`. Asserting that makes the day it starts working visible: this
    # test fails, and whatever adds pickling has to keep interning intact (see
    # `test_construction_returns_the_interned_member`) so that an unpickled
    # member is still the very object bound to the class attribute.
    member = members(enum_class)[0]
    with pytest.raises(TypeError):
        pickle.dumps(member)
