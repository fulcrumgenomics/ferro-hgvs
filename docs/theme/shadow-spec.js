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
 *   5. each executable row gets an expand toggle that reveals a biopython-style
 *      3-line alignment, fetched from the blessed alignments.json sidecar
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
            h2.classList.add("ss-clause"); // divider + clause-chip styling
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

    // --- Per-example alignments -------------------------------------------
    // Fetch the blessed sidecar (served at <shadow-spec root>/alignments.json), and for each
    // executable table row whose input has an entry, add a toggle that reveals a biopython-style
    // 3-line alignment in a row beneath it. Data-driven: rows without an entry get no toggle.

    function shadowRoot() {
        var p = window.location.pathname;
        var i = p.indexOf("/shadow-spec/");
        return i === -1 ? null : p.slice(0, i + "/shadow-spec/".length);
    }

    function esc(s) {
        return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
    }

    // One sequence line with its edited columns [es,ee) wrapped for the background tint.
    function hlLine(label, seq, es, ee, tint) {
        var body = tint
            ? esc(seq.slice(0, es)) +
              '<span class="ss-align-edit">' +
              esc(seq.slice(es, ee)) +
              "</span>" +
              esc(seq.slice(ee))
            : esc(seq);
        return '<span class="ss-align-gutter">' + label + "</span>" + body;
    }

    function alignBlock(a) {
        var pre =
            hlLine("ref", a.ref, a.editStart, a.editEnd, true) +
            "\n" +
            hlLine("   ", a.connector, a.editStart, a.editEnd, false) +
            "\n" +
            hlLine("alt", a.alt, a.editStart, a.editEnd, true);
        return (
            '<div class="ss-align">' +
            '<div class="ss-align-note">' +
            esc(a.note) +
            "</div>" +
            '<pre class="ss-align-pre">' +
            pre +
            "</pre>" +
            "</div>"
        );
    }

    var alignSeq = 0;

    function wireAlignments() {
        var root = shadowRoot();
        if (!root) return;
        fetch(root + "alignments.json")
            .then(function (r) {
                return r.ok ? r.json() : null;
            })
            .then(function (map) {
                if (!map) return;
                var tables = document.querySelectorAll(".content table.ss-examples");
                Array.prototype.forEach.call(tables, function (table) {
                    var nCols = table.querySelectorAll("thead th").length;
                    var rows = table.querySelectorAll("tbody tr");
                    Array.prototype.forEach.call(rows, function (tr) {
                        var inputCell = tr.querySelector(".ss-col-input");
                        var normCell = tr.querySelector(".ss-col-norm");
                        if (!inputCell || !normCell) return;
                        var key = inputCell.textContent.trim();
                        var a = map[key];
                        if (!a) return;

                        var id = "ss-aln-" + ++alignSeq;
                        var btn = document.createElement("button");
                        btn.className = "ss-align-toggle";
                        btn.type = "button";
                        btn.setAttribute("aria-expanded", "false");
                        btn.setAttribute("aria-controls", id);
                        var showLabel = "Show alignment for " + key;
                        var hideLabel = "Hide alignment for " + key;
                        btn.title = showLabel;
                        btn.setAttribute("aria-label", showLabel);
                        btn.innerHTML =
                            'align <span class="ss-align-caret" aria-hidden="true">▸</span>';

                        var alnRow = document.createElement("tr");
                        alnRow.className = "ss-align-row";
                        alnRow.id = id;
                        alnRow.hidden = true;
                        var td = document.createElement("td");
                        td.colSpan = nCols;
                        td.innerHTML = alignBlock(a);
                        alnRow.appendChild(td);
                        tr.parentNode.insertBefore(alnRow, tr.nextSibling);

                        btn.addEventListener("click", function () {
                            var open = alnRow.hidden;
                            alnRow.hidden = !open;
                            btn.setAttribute("aria-expanded", open ? "true" : "false");
                            btn.classList.toggle("ss-open", open);
                            btn.title = open ? hideLabel : showLabel;
                            btn.setAttribute("aria-label", open ? hideLabel : showLabel);
                            var caret = btn.querySelector(".ss-align-caret");
                            if (caret) caret.textContent = open ? "▾" : "▸";
                        });
                        normCell.appendChild(btn);
                    });
                });
            })
            .catch(function () {
                /* no alignments served (or offline) — leave the tables as-is */
            });
    }

    // The spec's MkDocs `!!! note "…"` admonitions are quoted verbatim (the source bytes are
    // pinned by the harness), so mdBook renders them raw: a literal `!!! note "…"` line plus the
    // indented answer as a code block. Reformat the RENDERED DOM to read like the note the spec
    // intends — display-only, the source is untouched.
    function inlineMarkdown(text) {
        // The answer text comes from the committed, byte-pinned spec quote (trusted), and may
        // already carry the spec's own HTML — notably <code class="invalid">…</code>, which the
        // spec uses to mark an invalid example and the theme styles. Pass that through, and render
        // the remaining Markdown inline (`code` and **bold**). Do NOT escape `<`/`>`, or the spec's
        // HTML shows as literal tag text.
        return text
            .replace(/`([^`]+)`/g, "<code>$1</code>")
            .replace(/\*\*([^*]+)\*\*/g, "<strong>$1</strong>");
    }

    function prettifySpecNotes() {
        var quotes = document.querySelectorAll(".content blockquote");
        Array.prototype.forEach.call(quotes, function (bq) {
            var first = bq.querySelector("p");
            if (!first) return;
            // Match a leading `!!! <type> "<title>"` line (mdBook curls the quotes).
            var m = first.innerHTML.match(/^\s*!!!\s*\w+\s*["“”]([\s\S]*)["“”]\s*$/);
            if (!m) return;
            bq.classList.add("ss-note");
            first.classList.add("ss-note-title");
            first.innerHTML = m[1]; // keep the title's already-rendered <code> spans
            var blocks = bq.querySelectorAll("pre > code");
            Array.prototype.forEach.call(blocks, function (code) {
                var p = document.createElement("p");
                p.className = "ss-note-body";
                p.innerHTML = inlineMarkdown(code.textContent);
                code.parentNode.replaceWith(p);
            });
        });
    }

    function run() {
        enhanceTables();
        buildOnPageNav();
        wireAlignments();
        prettifySpecNotes();
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", run);
    } else {
        run();
    }
})();
