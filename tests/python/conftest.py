"""Shared fixtures for the Python binding tests."""

import re
from collections.abc import Callable
from pathlib import Path

import pytest

# Location of the committed type stub, relative to this file
# (<repo>/tests/python/conftest.py -> <repo>/python/ferro_hgvs/__init__.pyi).
_STUB_PATH = Path(__file__).resolve().parents[2] / "python" / "ferro_hgvs" / "__init__.pyi"


def _stub_enum_members(class_name: str) -> dict[str, int]:
    """Parse the declared members of a native enum class from the type stub.

    Reads the block following ``class <class_name>(_NativeEnum):`` up to the
    next top-level ``class`` declaration, returning the declared name->value
    mapping. Members are annotated with their own class and carry the integer
    value in a trailing comment (``Plus: Strand  # 0``), because the value is
    not assignable under the annotation — see issue #1245, which moved these
    off the ``IntEnum`` declaration the runtime never honoured.

    Docstring lines and free-standing comments are ignored, so a retired
    discriminant documented only in a comment is not treated as a member.
    """
    lines = _STUB_PATH.read_text().splitlines()
    header = re.compile(rf"^class {re.escape(class_name)}\(_NativeEnum\):")
    member = re.compile(rf"^    (\w+): {re.escape(class_name)}\s+# (\d+)$")
    members: dict[str, int] = {}
    in_block = False
    for line in lines:
        if header.match(line):
            in_block = True
            continue
        if in_block:
            if line.startswith("class "):
                break
            match = member.match(line)
            if match:
                members[match.group(1)] = int(match.group(2))
    return members


@pytest.fixture
def stub_enum_members() -> Callable[[str], dict[str, int]]:
    """Return the stub-declared ``name -> value`` members of a native enum.

    Lives here rather than in one test module because two suites need it: the
    long-standing ``ErrorType`` drift guard (issue #630) and the all-enum guard
    added for #1245.
    """
    return _stub_enum_members
