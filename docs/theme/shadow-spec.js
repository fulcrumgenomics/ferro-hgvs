/**
 * Shadow-spec readability enhancer.
 *
 * The interpretation pages carry `| Input | Verdict | Normalizes to | Notes |`
 * tables whose cells are executed by the `shadow_spec` test: the Verdict cell
 * must contain the bare word `recommended` / `conformant` / `refused` / `bug`,
 * and the "Normalizes to" cell the bare token `self` or `—`. So the badges and
 * chips CANNOT be authored as inline HTML in the markdown (that would change the
 * cell text the harness matches). Instead we style the *rendered* DOM here:
 *
 *   1. verdict cells  -> colored badge span
 *   2. normalizes-to  -> `self`/`—` chips
 *   3. every example table gets `.ss-examples` + per-column classes so the
 *      compact/balanced column widths in custom.css apply regardless of the
 *      3- or 4-column shape
 *   4. an "On this page" clause index is built from the `file.md:lines` headings
 *
 * Purely cosmetic; the source bytes the harness reads are untouched. Scoped to
 * /shadow-spec/ pages.
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
    var COL_CLASS = {
        input: "ss-col-input",
        spelling: "ss-col-input",
        verdict: "ss-col-verdict",
        "normalizes to": "ss-col-norm",
        notes: "ss-col-notes",
    };

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

    function enhanceTables() {
        var tables = document.querySelectorAll(".content table");
        Array.prototype.forEach.call(tables, function (table) {
            var headCells = table.querySelectorAll("thead th");
            var headers = Array.prototype.map.call(headCells, function (th) {
                return th.textContent.trim().toLowerCase();
            });
            var vIdx = headers.indexOf("verdict");
            if (vIdx === -1) {
                return; // not an example table
            }
            table.classList.add("ss-examples");
            var nIdx = headers.indexOf("normalizes to");

            // tag every column (header + body) so CSS widths hold for 3/4 cols
            var colClass = headers.map(function (h) {
                return COL_CLASS[h] || "";
            });
            Array.prototype.forEach.call(headCells, function (th, i) {
                if (colClass[i]) th.classList.add(colClass[i]);
            });
            var rows = table.querySelectorAll("tbody tr");
            Array.prototype.forEach.call(rows, function (tr) {
                var cells = tr.children;
                Array.prototype.forEach.call(cells, function (td, i) {
                    if (colClass[i]) td.classList.add(colClass[i]);
                });
                if (cells[vIdx]) {
                    var word = cells[vIdx].textContent.trim().toLowerCase();
                    var b = badge(word);
                    if (b) {
                        cells[vIdx].textContent = "";
                        cells[vIdx].appendChild(b);
                    }
                }
                if (nIdx !== -1 && cells[nIdx]) {
                    var chip = normChip(cells[nIdx].textContent.trim());
                    if (chip) {
                        cells[nIdx].textContent = "";
                        cells[nIdx].appendChild(chip);
                    }
                }
            });
        });
    }

    // A clause heading renders as <h2 id="..."><code>file.md:lines</code> — title</h2>.
    var CLAUSE_RE = /^[\w./-]+\.md:[0-9-]+$/;

    function buildOnPageNav() {
        var content = document.querySelector(".content main") || document.querySelector(".content");
        if (!content) return;
        var h2s = content.querySelectorAll("h2");
        var entries = [];
        Array.prototype.forEach.call(h2s, function (h2) {
            var code = h2.querySelector("code");
            if (!code || !CLAUSE_RE.test(code.textContent.trim()) || !h2.id) return;
            var token = code.textContent.trim();
            var line = token.slice(token.indexOf(":")); // ":5", ":16-17"
            // title = heading text after the em dash, trimmed
            var full = h2.textContent.trim();
            var title = full;
            var dash = full.indexOf("—");
            if (dash !== -1) title = full.slice(dash + 1).trim();
            entries.push({ id: h2.id, line: line, title: title, h2: h2 });
        });
        if (entries.length < 3) return; // not worth a nav on tiny pages

        var nav = document.createElement("nav");
        nav.className = "ss-onpage-nav";
        var t = document.createElement("span");
        t.className = "ss-nav-title";
        t.textContent = "On this page";
        nav.appendChild(t);
        entries.forEach(function (e) {
            var a = document.createElement("a");
            a.href = "#" + e.id;
            var ln = document.createElement("span");
            ln.className = "ss-nav-line";
            ln.textContent = e.line;
            a.appendChild(ln);
            a.appendChild(document.createTextNode(" " + e.title));
            nav.appendChild(a);
        });
        // insert just before the first clause heading
        var first = entries[0].h2;
        first.parentNode.insertBefore(nav, first);
    }

    function run() {
        enhanceTables();
        buildOnPageNav();
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", run);
    } else {
        run();
    }
})();
