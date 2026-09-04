/**
 * Shadow-spec readability enhancer.
 *
 * The interpretation pages carry `| Input | Verdict | Normalizes to | Notes |`
 * tables whose cells are executed by the `shadow_spec` test: the Verdict cell
 * must contain the bare word `recommended` / `conformant` / `refused` / `bug`,
 * and the "Normalizes to" cell the bare token `self` or `—`. So the badges and
 * chips CANNOT be authored as inline HTML in the markdown (that would change the
 * cell text the harness matches). Instead we style the *rendered* DOM here:
 * find each table's Verdict / Normalizes-to columns and wrap the token in a
 * styled span. Purely cosmetic; the source bytes are untouched.
 *
 * Scoped to /shadow-spec/ pages so it never touches tables elsewhere.
 */
(function () {
    if (!/\/shadow-spec\//.test(window.location.pathname)) {
        return;
    }

    var VERDICT_CLASS = {
        recommended: "verdict-recommended",
        conformant: "verdict-conformant",
        refused: "verdict-refused",
        bug: "verdict-bug",
    };

    function headerIndex(table, name) {
        var ths = table.querySelectorAll("thead th");
        for (var i = 0; i < ths.length; i++) {
            if (ths[i].textContent.trim().toLowerCase() === name) {
                return i;
            }
        }
        return -1;
    }

    function badge(word) {
        var cls = VERDICT_CLASS[word];
        if (!cls) {
            return null;
        }
        var span = document.createElement("span");
        span.className = "verdict-badge " + cls;
        span.textContent = word;
        return span;
    }

    function normChip(raw) {
        var span = document.createElement("span");
        if (raw === "self") {
            span.className = "norm-chip norm-self";
            span.textContent = "= input";
            span.title =
                "Fixed point: ferro leaves this input unchanged (normalize(input) == input).";
            return span;
        }
        if (raw === "—" || raw === "-") {
            span.className = "norm-chip norm-na";
            span.textContent = "not run here";
            span.title =
                "Not run against the test reference on this page (usually an accession outside " +
                "the committed slice). The row is still parse-checked.";
            return span;
        }
        return null;
    }

    function enhance() {
        var tables = document.querySelectorAll(".content table");
        Array.prototype.forEach.call(tables, function (table) {
            var vIdx = headerIndex(table, "verdict");
            if (vIdx === -1) {
                return;
            }
            var nIdx = headerIndex(table, "normalizes to");
            var rows = table.querySelectorAll("tbody tr");
            Array.prototype.forEach.call(rows, function (tr) {
                var cells = tr.children;
                if (cells[vIdx]) {
                    var word = cells[vIdx].textContent.trim().toLowerCase();
                    var b = badge(word);
                    if (b) {
                        cells[vIdx].textContent = "";
                        cells[vIdx].appendChild(b);
                    }
                }
                if (nIdx !== -1 && cells[nIdx]) {
                    var raw = cells[nIdx].textContent.trim();
                    var chip = normChip(raw);
                    if (chip) {
                        cells[nIdx].textContent = "";
                        cells[nIdx].appendChild(chip);
                    }
                }
            });
        });
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", enhance);
    } else {
        enhance();
    }
})();
