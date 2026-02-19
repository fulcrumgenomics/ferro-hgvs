// HGVS Normalizer Web UI

let availableTools = [];

// Help system with 6 tabs
class HelpSystem {
    constructor() {
        this.content = {
            overview: `
                <h4>What is ferro-hgvs?</h4>
                <p>A multi-tool web service for normalizing HGVS variant nomenclature.
                Compare results from up to four different normalization tools to validate variants
                and identify discrepancies.</p>
                <h4>Quick Start</h4>
                <ol>
                    <li>Select which tools to use (checkboxes above)</li>
                    <li>Choose an error mode (Silent, Lenient, or Strict)</li>
                    <li>Enter an HGVS variant expression</li>
                    <li>Click "Normalize" to see results from all selected tools</li>
                </ol>
                <h4>Features</h4>
                <ul>
                    <li>Single variant and batch processing modes</li>
                    <li>Side-by-side comparison of multiple normalization tools</li>
                    <li>Detailed component breakdown for validated variants</li>
                    <li>Agreement detection to identify tool discrepancies</li>
                </ul>
            `,
            tools: `
                <h4>Available Tools</h4>
                <dl>
                    <dt><strong>ferro</strong> <a href="https://github.com/fulcrumgenomics/ferro-hgvs" target="_blank" rel="noopener">[GitHub]</a></dt>
                    <dd>Native Rust implementation of HGVS parsing and normalization. Fast, local-only processing
                    with no external dependencies. Most comprehensive normalization support including intronic variants (c.117-2del) that other tools cannot normalize.</dd>

                    <dt><strong>mutalyzer</strong> <a href="https://mutalyzer.nl/" target="_blank" rel="noopener">[Website]</a> <a href="https://github.com/mutalyzer/mutalyzer" target="_blank" rel="noopener">[GitHub]</a></dt>
                    <dd>Uses the Mutalyzer service (via API or local Python subprocess). Supports extensive
                    variant types including complex genomic rearrangements. Requires network access or local Python installation.</dd>

                    <dt><strong>biocommons</strong> <a href="https://github.com/biocommons/hgvs" target="_blank" rel="noopener">[GitHub]</a> <a href="https://hgvs.readthedocs.io/" target="_blank" rel="noopener">[Docs]</a></dt>
                    <dd>Python biocommons/hgvs library. The reference implementation for HGVS normalization.
                    Requires UTA database and SeqRepo for full functionality.</dd>

                    <dt><strong>hgvs-rs</strong> <a href="https://github.com/varfish-org/hgvs-rs" target="_blank" rel="noopener">[GitHub]</a></dt>
                    <dd>Native Rust port of the biocommons library. Provides similar functionality with
                    better performance. Same data requirements as biocommons (UTA + SeqRepo).</dd>
                </dl>
                <h4>Tool Selection</h4>
                <p>Tools that are unavailable (grayed out) are not configured or their dependencies are not accessible.
                Check the service status in the footer for more details.</p>
            `,
            support: `
                <h4>HGVS Syntax &amp; Tool Support</h4>
                <p>HGVS (<a href="https://hgvs-nomenclature.org/" target="_blank" rel="noopener">Human Genome Variation Society</a>) nomenclature
                describes sequence variants. Format: <code>Reference:Coordinate.Change</code> (e.g., <code>NM_000088.3:c.589G>T</code>)</p>
                <p>This matrix shows which features each tool supports: <strong>V/N</strong> = Validate &amp; Normalize,
                <strong>V</strong> = Validate only, <strong>N</strong> = Normalize only, <strong>-</strong> = Not supported.</p>

                <h5>Reference Types <a href="https://hgvs-nomenclature.org/stable/background/refseq/" target="_blank" rel="noopener">[Spec]</a></h5>
                <table>
                    <tr><th>Reference</th><th>Description</th><th>Example</th><th>ferro</th><th>mutalyzer</th><th>biocommons</th><th>hgvs-rs</th></tr>
                    <tr><td><code>NM_</code></td><td>Coding transcript (mRNA)</td><td><code>NM_000088.3</code></td><td>V/N</td><td>V/N</td><td>V/N</td><td>V/N</td></tr>
                    <tr><td><code>NR_</code></td><td>Non-coding transcript</td><td><code>NR_024540.1</code></td><td>V/N</td><td>V/N</td><td>V/N</td><td>V/N</td></tr>
                    <tr><td><code>NC_</code></td><td>Genomic chromosome</td><td><code>NC_000001.11</code></td><td>V/N</td><td>V/N</td><td>V/N</td><td>V/N</td></tr>
                    <tr><td><code>NG_</code></td><td>Genomic gene region</td><td><code>NG_007400.1</code></td><td>V/N</td><td>V/N</td><td>V/N</td><td>V/N</td></tr>
                    <tr><td><code>NP_</code></td><td>Protein</td><td><code>NP_000079.2</code></td><td>V</td><td>V/N</td><td>V/N</td><td>-</td></tr>
                    <tr><td><code>LRG_</code></td><td>Locus Reference Genomic</td><td><code>LRG_1</code></td><td>V/N</td><td>V/N</td><td>V/N</td><td>V/N</td></tr>
                    <tr><td><code>ENST</code></td><td>Ensembl transcript</td><td><code>ENST00000357033</code></td><td>-</td><td>V/N</td><td>V/N</td><td>-</td></tr>
                </table>

                <h5>Coordinate Types</h5>
                <table>
                    <tr><th>Coordinate</th><th>Description</th><th>Example</th><th>ferro</th><th>mutalyzer</th><th>biocommons</th><th>hgvs-rs</th></tr>
                    <tr><td><code>c.</code></td><td><a href="https://hgvs-nomenclature.org/stable/recommendations/DNA/variant/substitution/" target="_blank" rel="noopener">Coding DNA</a> (relative to CDS)</td><td><code>c.589G>T</code></td><td>V/N</td><td>V/N</td><td>V/N</td><td>V/N</td></tr>
                    <tr><td><code>g.</code></td><td><a href="https://hgvs-nomenclature.org/stable/recommendations/DNA/variant/substitution/" target="_blank" rel="noopener">Genomic</a> (absolute position)</td><td><code>g.12345A>G</code></td><td>V/N</td><td>V/N</td><td>V/N</td><td>V/N</td></tr>
                    <tr><td><code>n.</code></td><td><a href="https://hgvs-nomenclature.org/stable/recommendations/DNA/variant/substitution/" target="_blank" rel="noopener">Non-coding transcript</a></td><td><code>n.100del</code></td><td>V/N</td><td>V/N</td><td>V/N</td><td>V/N</td></tr>
                    <tr><td><code>p.</code></td><td><a href="https://hgvs-nomenclature.org/stable/recommendations/protein/" target="_blank" rel="noopener">Protein</a></td><td><code>p.Gly12Val</code></td><td>V</td><td>V/N</td><td>V/N</td><td>-</td></tr>
                    <tr><td><code>m.</code></td><td><a href="https://hgvs-nomenclature.org/stable/recommendations/DNA/variant/substitution/" target="_blank" rel="noopener">Mitochondrial</a></td><td><code>m.8993T>G</code></td><td>V/N</td><td>V/N</td><td>V/N</td><td>V/N</td></tr>
                    <tr><td><code>r.</code></td><td><a href="https://hgvs-nomenclature.org/stable/recommendations/RNA/" target="_blank" rel="noopener">RNA</a></td><td><code>r.76a>u</code></td><td>V</td><td>V/N</td><td>V/N</td><td>-</td></tr>
                    <tr><td><code>+/-</code></td><td>Intronic positions</td><td><code>c.100+5G>A</code></td><td>V/N</td><td>V/N</td><td>V/N</td><td>V</td></tr>
                </table>

                <h5>Variant Types <a href="https://hgvs-nomenclature.org/stable/recommendations/DNA/variant/" target="_blank" rel="noopener">[Spec]</a></h5>
                <table>
                    <tr><th>Type</th><th>Description</th><th>Example</th><th>ferro</th><th>mutalyzer</th><th>biocommons</th><th>hgvs-rs</th></tr>
                    <tr><td><a href="https://hgvs-nomenclature.org/stable/recommendations/DNA/variant/substitution/" target="_blank" rel="noopener">Substitution</a></td><td>Single nucleotide change</td><td><code>c.589G>T</code></td><td>V/N</td><td>V/N</td><td>V/N</td><td>V/N</td></tr>
                    <tr><td><a href="https://hgvs-nomenclature.org/stable/recommendations/DNA/variant/deletion/" target="_blank" rel="noopener">Deletion</a></td><td>Nucleotide(s) removed</td><td><code>c.589del</code></td><td>V/N</td><td>V/N</td><td>V/N</td><td>V/N</td></tr>
                    <tr><td><a href="https://hgvs-nomenclature.org/stable/recommendations/DNA/variant/insertion/" target="_blank" rel="noopener">Insertion</a></td><td>Nucleotide(s) added</td><td><code>c.589_590insA</code></td><td>V/N</td><td>V/N</td><td>V/N</td><td>V/N</td></tr>
                    <tr><td><a href="https://hgvs-nomenclature.org/stable/recommendations/DNA/variant/duplication/" target="_blank" rel="noopener">Duplication</a></td><td>Sequence copied in tandem</td><td><code>c.589dup</code></td><td>V/N</td><td>V/N</td><td>V/N</td><td>V/N</td></tr>
                    <tr><td><a href="https://hgvs-nomenclature.org/stable/recommendations/DNA/variant/delins/" target="_blank" rel="noopener">Delins</a></td><td>Deletion + insertion</td><td><code>c.589delinsAT</code></td><td>V/N</td><td>V/N</td><td>V/N</td><td>V/N</td></tr>
                    <tr><td><a href="https://hgvs-nomenclature.org/stable/recommendations/DNA/variant/inversion/" target="_blank" rel="noopener">Inversion</a></td><td>Sequence reversed</td><td><code>c.589_600inv</code></td><td>V</td><td>V/N</td><td>V/N</td><td>V</td></tr>
                    <tr><td><a href="https://hgvs-nomenclature.org/stable/recommendations/DNA/variant/repeated/" target="_blank" rel="noopener">Repeat</a></td><td>Tandem repeat expansion</td><td><code>c.589CAG[23]</code></td><td>V</td><td>V/N</td><td>V</td><td>-</td></tr>
                    <tr><td>Conversion</td><td>Gene conversion event</td><td><code>c.589_600con</code></td><td>-</td><td>V/N</td><td>V</td><td>-</td></tr>
                </table>

                <p><em>Note: Actual support may vary by tool version and configuration. "Validate" = parsing/syntax checking,
                "Normalize" = 3' shifting and HGVS-compliant formatting.</em></p>
            `,
            operations: `
                <h4>Available Operations</h4>
                <dl>
                    <dt><strong>Normalize</strong></dt>
                    <dd>Validate and normalize HGVS variants using multiple tools. Compare results across ferro, mutalyzer, biocommons, and hgvs-rs to identify discrepancies.</dd>

                    <dt><strong>Convert</strong></dt>
                    <dd>Convert HGVS variants between coordinate systems:
                        <ul>
                            <li><code>c.</code> (coding) ↔ <code>g.</code> (genomic) - requires transcript data</li>
                            <li><code>c.</code> (coding) → <code>p.</code> (protein) - calculates amino acid position</li>
                            <li><code>n.</code> (non-coding) ↔ <code>g.</code> (genomic)</li>
                        </ul>
                        Requires cdot transcript data to be configured.
                    </dd>

                    <dt><strong>Effect</strong></dt>
                    <dd>Predict variant consequences using Sequence Ontology (SO) terms:
                        <ul>
                            <li>Classifies variants as splice_site_variant, frameshift_variant, inframe_deletion, etc.</li>
                            <li>Reports impact level: HIGH, MODERATE, LOW, MODIFIER</li>
                            <li>Optional NMD (Nonsense-Mediated Decay) prediction for truncating variants</li>
                            <li>Protein consequence with amino acid position</li>
                        </ul>
                    </dd>

                    <dt><strong>Liftover</strong></dt>
                    <dd>Convert genomic coordinates between genome builds:
                        <ul>
                            <li>GRCh37 (hg19) → GRCh38 (hg38)</li>
                            <li>GRCh38 (hg38) → GRCh37 (hg19)</li>
                        </ul>
                        Accepts positions as <code>chr7:117120148</code> or HGVS <code>NC_000007.13:g.117120148</code>.
                        Requires liftover chain files to be configured.
                    </dd>

                    <dt><strong>VCF</strong></dt>
                    <dd>Bidirectional conversion between VCF and HGVS formats:
                        <ul>
                            <li><strong>HGVS → VCF:</strong> Convert genomic HGVS (g.) to VCF CHROM/POS/REF/ALT</li>
                            <li><strong>VCF → HGVS:</strong> Convert VCF fields to HGVS g. notation, optionally with c./p. if transcript provided</li>
                        </ul>
                        Note: Indels use "N" placeholder for padding bases when sequence data unavailable.
                    </dd>
                </dl>
            `,
            api: `
                <h4>REST API</h4>
                <h5>Normalize Single Variant</h5>
<pre><code>POST /api/v1/normalize
Content-Type: application/json

{
  "hgvs": "NM_000088.3:c.589G>T",
  "tools": ["ferro", "mutalyzer"],
  "error_mode": "lenient"
}</code></pre>
                <h5>Batch Normalize</h5>
<pre><code>POST /api/v1/batch/normalize
Content-Type: application/json

{
  "variants": ["NM_000088.3:c.589G>T", "NM_000088.3:c.590del"],
  "tools": ["ferro"],
  "error_mode": "strict"
}</code></pre>
                <h5>Validate Variant (ferro only)</h5>
<pre><code>POST /api/v1/validate
Content-Type: application/json

{
  "hgvs": "NM_000088.3:c.589G>T"
}</code></pre>
                <h5>Health Check</h5>
<pre><code>GET /health</code></pre>

                <h5>Convert Coordinates</h5>
<pre><code>POST /api/v1/convert
Content-Type: application/json

{
  "hgvs": "NM_000249.4:c.350C>T",
  "target_system": "g",
  "include_all": false
}</code></pre>

                <h5>Effect Prediction</h5>
<pre><code>POST /api/v1/effect
Content-Type: application/json

{
  "hgvs": "NM_000249.4:c.350C>T",
  "include_nmd": true
}</code></pre>

                <h5>Liftover</h5>
<pre><code>POST /api/v1/liftover
Content-Type: application/json

{
  "position": "chr7:117120148",
  "from_build": "GRCh37",
  "to_build": "GRCh38"
}</code></pre>

                <h5>HGVS to VCF</h5>
<pre><code>POST /api/v1/hgvs-to-vcf
Content-Type: application/json

{
  "hgvs": "NC_000007.14:g.117559593G>A",
  "build": "GRCh38"
}</code></pre>

                <h5>VCF to HGVS</h5>
<pre><code>POST /api/v1/vcf-to-hgvs
Content-Type: application/json

{
  "chrom": "chr7",
  "pos": 117559593,
  "ref": "G",
  "alt": "A",
  "build": "GRCh38",
  "transcript": "NM_000249.4"  // optional
}</code></pre>
                <h4>Response Format</h4>
                <p>All endpoints return JSON with results, processing time, and any errors or warnings.</p>
            `,
            errors: `
                <h4>Error Modes</h4>
                <p>Error modes control how the service handles errors during normalization.</p>
                <dl>
                    <dt><strong>Silent</strong></dt>
                    <dd>Suppress all errors and warnings. Returns empty or null results on failure.
                    Useful when you only care about successful normalizations and want to ignore failures.</dd>

                    <dt><strong>Lenient</strong> (default)</dt>
                    <dd>Log warnings but continue processing. Returns partial results when possible.
                    Best for batch operations where you want to process as many variants as possible.</dd>

                    <dt><strong>Strict</strong></dt>
                    <dd>Fail immediately on any error. Returns detailed error information.
                    Best for validation workflows where you need to know about every issue.</dd>
                </dl>
                <h4>Error Categories</h4>
                <ul>
                    <li><strong>Parse errors</strong> - Invalid HGVS syntax or format</li>
                    <li><strong>Reference errors</strong> - Sequence or transcript not found</li>
                    <li><strong>Validation errors</strong> - Position out of bounds, unsupported variant type</li>
                    <li><strong>Tool errors</strong> - Tool unavailable or execution failed</li>
                </ul>
            `,
            examples: `
                <h4>Example Variants (work with all 4 tools)</h4>
                <table>
                    <tr><th>Type</th><th>Example</th><th>Description</th></tr>
                    <tr><td>Substitution</td><td><code>NM_000249.4:c.350C>T</code></td><td>C to T at position 350 (MLH1)</td></tr>
                    <tr><td>Deletion</td><td><code>NM_000249.4:c.1852_1853delAA</code></td><td>Delete AA at 1852-1853 (MLH1)</td></tr>
                    <tr><td>Deletion</td><td><code>NM_007294.4:c.68_69delAG</code></td><td>BRCA1 185delAG founder mutation</td></tr>
                    <tr><td>Duplication</td><td><code>NM_007294.4:c.5266dupC</code></td><td>BRCA1 5382insC founder mutation</td></tr>
                </table>
                <h4>Intronic Variants (ferro only)</h4>
                <table>
                    <tr><th>Type</th><th>Example</th><th>Description</th></tr>
                    <tr><td>Intronic</td><td><code>NM_000249.4:c.117-2del</code></td><td>Splice acceptor deletion (MLH1)</td></tr>
                    <tr><td>Intronic</td><td><code>NM_000249.4:c.116+1G>A</code></td><td>Splice donor substitution (MLH1)</td></tr>
                </table>
                <p><em>Note: Only ferro supports intronic variant normalization. Other tools require genomic coordinates.</em></p>
                <h4>Normalization Example</h4>
                <p>Input: <code>NM_000249.4:c.1852_1853delAA</code></p>
                <p>Normalized: <code>NM_000249.4:c.1852_1853del</code></p>
                <p>The explicit deleted sequence (AA) is removed per HGVS guidelines.</p>
            `
        };
    }

    init() {
        const toggle = document.querySelector('.help-toggle');
        const content = document.querySelector('.help-content');
        const tabs = document.querySelectorAll('.help-tab-btn');

        if (toggle && content) {
            toggle.addEventListener('click', () => {
                const expanded = toggle.getAttribute('aria-expanded') === 'true';
                toggle.setAttribute('aria-expanded', !expanded);
                content.hidden = expanded;
                if (!expanded) {
                    this.showTab('overview');
                }
            });
        }

        tabs.forEach(tab => {
            tab.addEventListener('click', () => this.showTab(tab.dataset.tab));
        });
    }

    showTab(tabId) {
        document.querySelectorAll('.help-tab-btn').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.tab === tabId);
        });
        const container = document.getElementById('help-tab-content');
        if (container && this.content[tabId]) {
            container.innerHTML = this.content[tabId];
        }
    }
}

// Initialize on page load
document.addEventListener('DOMContentLoaded', () => {
    new HelpSystem().init();
    refreshHealth();

    // Poll for updated health cache every 60 seconds
    // (Server runs actual health checks every 15 minutes in background)
    setInterval(refreshHealth, 60000);
});

// Refresh health status (quick check)
async function refreshHealth() {
    const indicator = document.getElementById('status-indicator');
    const text = document.getElementById('status-text');
    const healthSummary = document.getElementById('health-summary');

    // Show loading state only on first load (when indicator has no class)
    const isFirstLoad = indicator && !indicator.classList.contains('healthy') &&
                        !indicator.classList.contains('unhealthy') &&
                        !indicator.classList.contains('partial');

    if (isFirstLoad) {
        if (indicator) indicator.className = 'status-indicator';
        if (text) text.textContent = 'Loading health status...';
        if (healthSummary) healthSummary.innerHTML = '<div class="health-loading">Loading health status...</div>';
    }

    try {
        const response = await fetch('/health');
        const data = await response.json();

        availableTools = data.available_tools || [];

        // Update footer with basic status
        updateFooterStatus(data);

        // Update tool checkboxes
        updateToolCheckboxes(data.tools || []);

        // Fetch detailed health (reads cached results from server)
        fetchDetailedHealth();

    } catch (error) {
        if (indicator) indicator.className = 'status-indicator unhealthy';
        if (text) text.textContent = 'Service unavailable';
    }
}

// Fetch detailed health check results (from server cache)
async function fetchDetailedHealth() {
    const container = document.getElementById('health-summary');

    try {
        const response = await fetch('/health/detailed');
        const data = await response.json();
        const testResults = data.test_results || [];

        if (testResults.length === 0) {
            // Cache not populated yet - server is still running initial health check
            if (container) {
                container.innerHTML = '<div class="health-loading">Health checks running... (first startup)</div>';
            }
        } else {
            updateHealthSummary(testResults);
        }
    } catch (error) {
        console.error('Failed to fetch detailed health:', error);
        if (container) {
            container.innerHTML = '<div class="health-error">Failed to load health check</div>';
        }
    }
}

// Update footer with basic status
function updateFooterStatus(healthData) {
    const indicator = document.getElementById('status-indicator');
    const text = document.getElementById('status-text');

    if (!indicator || !text) return;

    const availableCount = healthData.available_tools?.length || 0;

    if (healthData.status === 'starting') {
        indicator.className = 'status-indicator';
        text.textContent = 'Starting up...';
    } else if (healthData.status === 'healthy') {
        indicator.className = 'status-indicator healthy';
        text.textContent = 'Service healthy';
    } else if (availableCount > 0) {
        indicator.className = 'status-indicator partial';
        text.textContent = 'Partial availability';
    } else {
        indicator.className = 'status-indicator unhealthy';
        text.textContent = 'Service unavailable';
    }
}

// Update health summary with test results
function updateHealthSummary(testResults) {
    const container = document.getElementById('health-summary');
    if (!container) return;

    // Fixed tool order - always show all 4
    const toolOrder = ['ferro', 'mutalyzer', 'biocommons', 'hgvs-rs'];

    // Build a map of results by tool name
    const resultsByTool = {};
    for (const result of testResults) {
        resultsByTool[result.tool] = result;
    }

    // Get max total_tests (all patterns) for color gradient baseline
    const maxTotalTests = testResults.reduce((max, r) => Math.max(max, r.total_tests || 0), 0);

    // Helper: calculate color gradient (green to red based on percentage)
    function getGradientColor(value, max) {
        if (max === 0) return '#6c757d';  // Grey for no data
        const percent = value / max;
        // Green (120) to Red (0) in HSL
        const hue = Math.round(percent * 120);
        return `hsl(${hue}, 70%, 45%)`;
    }

    // Build summary row in fixed order
    let summaryHtml = '<div class="health-summary-row" onclick="toggleHealthDetails()">';

    for (const toolName of toolOrder) {
        const result = resultsByTool[toolName];
        if (result && result.total > 0) {
            // Tool has applicable tests
            const healthPercent = Math.round((result.passed / result.total) * 100);
            const coverageCount = result.total;

            // Status based on health rate (of things it supports, are they healthy?)
            const statusClass = healthPercent === 100 ? 'pass' : (healthPercent >= 50 ? 'partial' : 'fail');
            const icon = healthPercent === 100 ? '✓' : (healthPercent >= 50 ? '⚠' : '✗');
            const modeStr = result.mode ? ` (${result.mode})` : '';

            // Color gradient based on coverage (how many of max patterns does this tool support?)
            const coverageColor = getGradientColor(coverageCount, maxTotalTests);

            // Hybrid display: health as percentage, coverage as fraction
            summaryHtml += `<span class="tool-health ${statusClass}">
                <span class="tool-health-icon">${icon}</span>
                <span class="tool-health-name">${toolName}${modeStr}</span>
                <span class="tool-health-metrics">
                    <span class="health-rate" title="Health: ${result.passed}/${result.total} supported tests pass">${healthPercent}%</span>
                    <span class="coverage-rate" style="color: ${coverageColor}" title="Coverage: ${coverageCount}/${maxTotalTests} total patterns supported">(${coverageCount}/${maxTotalTests})</span>
                </span>
            </span>`;
        } else if (result && result.total === 0) {
            // Tool exists but all tests are N/A - grey out
            const modeStr = result.mode ? ` (${result.mode})` : '';
            summaryHtml += `<span class="tool-health unavailable">
                <span class="tool-health-icon">–</span>
                <span class="tool-health-name">${toolName}${modeStr}</span>
                <span class="tool-health-metrics">
                    <span class="health-rate">0%</span>
                    <span class="coverage-rate" title="Coverage: 0/${maxTotalTests} patterns supported">(0/${maxTotalTests})</span>
                </span>
            </span>`;
        } else {
            // Tool not available at all - grey out
            summaryHtml += `<span class="tool-health unavailable">
                <span class="tool-health-icon">–</span>
                <span class="tool-health-name">${toolName}</span>
                <span class="tool-health-metrics">N/A</span>
            </span>`;
        }
    }

    summaryHtml += '<span class="expand-icon">▼</span></div>';

    // Build unified details table with tools as columns (hidden by default)
    summaryHtml += '<div class="health-details" id="health-details" style="display: none;">';

    // Collect all unique tests across all tools
    const allTests = [];
    const testMap = {}; // Map of "category|name|variant" -> {category, name, variant, results: {tool: result}}

    for (const toolName of toolOrder) {
        const result = resultsByTool[toolName];
        if (!result) continue;

        for (const category of result.categories) {
            for (const test of category.tests) {
                const key = `${category.name}|${test.name}|${test.variant}`;
                if (!testMap[key]) {
                    testMap[key] = {
                        category: category.name,
                        name: test.name,
                        variant: test.variant,
                        results: {}
                    };
                    allTests.push(testMap[key]);
                }
                testMap[key].results[toolName] = test;
            }
        }
    }

    // Build table with tools as columns
    summaryHtml += '<table class="health-matrix-table">';

    // Header row with tool names
    summaryHtml += '<thead><tr>';
    summaryHtml += '<th>Category</th><th>Test</th><th>Variant</th>';
    for (const toolName of toolOrder) {
        const result = resultsByTool[toolName];
        const modeStr = result?.mode ? ` (${result.mode})` : '';
        const available = !!result;
        summaryHtml += `<th class="${available ? '' : 'unavailable'}">${toolName}${modeStr}</th>`;
    }
    summaryHtml += '</tr></thead>';

    // Group tests by category
    let currentCategory = '';
    summaryHtml += '<tbody>';

    for (const test of allTests) {
        const showCategory = test.category !== currentCategory;
        currentCategory = test.category;

        summaryHtml += '<tr>';
        summaryHtml += `<td class="category-cell">${showCategory ? escapeHtml(test.category) : ''}</td>`;
        summaryHtml += `<td class="test-name-cell">${escapeHtml(test.name)}</td>`;
        summaryHtml += `<td class="variant-cell"><code>${escapeHtml(test.variant)}</code></td>`;

        for (const toolName of toolOrder) {
            const toolResult = test.results[toolName];
            if (!resultsByTool[toolName]) {
                // Tool not available
                summaryHtml += '<td class="result-cell unavailable">–</td>';
            } else if (toolResult) {
                // Use status field: pass, fail, or na
                const status = toolResult.status || (toolResult.passed ? 'pass' : 'fail');
                let icon, statusClass;
                if (status === 'pass') {
                    icon = '✓';
                    statusClass = 'pass';
                } else if (status === 'na') {
                    icon = '⚠';
                    statusClass = 'na';
                } else {
                    icon = '✗';
                    statusClass = 'fail';
                }
                const tooltip = toolResult.error ? ` title="${escapeHtml(toolResult.error)}"` : '';
                summaryHtml += `<td class="result-cell ${statusClass}"${tooltip}>${icon}</td>`;
            } else {
                summaryHtml += '<td class="result-cell">–</td>';
            }
        }

        summaryHtml += '</tr>';
    }

    summaryHtml += '</tbody></table>';
    summaryHtml += '</div>';

    container.innerHTML = summaryHtml;
}

// Toggle health details visibility
function toggleHealthDetails() {
    const details = document.getElementById('health-details');
    const icon = document.querySelector('.expand-icon');
    if (details) {
        const isHidden = details.style.display === 'none';
        details.style.display = isHidden ? 'block' : 'none';
        if (icon) icon.textContent = isHidden ? '▲' : '▼';
    }
}

// Update tool selection checkboxes to show availability and mode
function updateToolCheckboxes(tools) {
    const allToolNames = ['ferro', 'mutalyzer', 'biocommons', 'hgvs-rs'];

    // Build maps of tool availability and mode from health check
    const toolAvailability = {};
    const toolMode = {};
    for (const tool of tools) {
        toolAvailability[tool.tool] = tool.available;
        toolMode[tool.tool] = tool.mode;
    }

    // Update each checkbox
    for (const toolName of allToolNames) {
        const container = document.querySelector(`[data-tool="${toolName}"]`);
        if (!container) continue;

        const checkbox = container.querySelector('input[type="checkbox"]');
        const nameSpan = container.querySelector('.tool-name');
        const available = toolAvailability[toolName] === true;
        const mode = toolMode[toolName];

        container.classList.toggle('unavailable', !available);

        if (checkbox) {
            checkbox.disabled = !available;
            if (!available) {
                checkbox.checked = false;
            }
        }

        // Show mode indicator for tools that have modes (e.g., mutalyzer)
        if (nameSpan && mode) {
            // Remove any existing mode indicator
            const existingMode = container.querySelector('.tool-mode');
            if (existingMode) existingMode.remove();

            // Add mode indicator
            const modeSpan = document.createElement('span');
            modeSpan.className = 'tool-mode';
            modeSpan.textContent = `(${mode})`;
            modeSpan.title = mode === 'api' ? 'Using remote API' : 'Using local Python subprocess';
            nameSpan.parentNode.insertBefore(modeSpan, nameSpan.nextSibling);
        }
    }
}

// Get selected tools (only enabled/available ones)
function getSelectedTools() {
    const checkboxes = document.querySelectorAll('#tool-checkboxes input[type="checkbox"]:checked:not(:disabled)');
    return Array.from(checkboxes).map(cb => cb.value);
}

// Get selected error mode
function getErrorMode() {
    const selected = document.querySelector('input[name="error_mode"]:checked');
    return selected ? selected.value : 'lenient';
}

// Show input mode (single/batch)
function showInputMode(mode) {
    document.querySelectorAll('.input-mode').forEach(el => el.classList.remove('active'));
    document.querySelectorAll('.input-mode-tabs .tab-button').forEach(el => el.classList.remove('active'));

    document.getElementById(`mode-${mode}`).classList.add('active');
    event.target.classList.add('active');
}

// Set example variant
function setExample(variant) {
    document.getElementById('hgvs-input').value = variant;
}

// Load batch example
function loadBatchExample() {
    document.getElementById('batch-input').value =
        'NM_000249.4:c.350C>T\nNM_000249.4:c.1852_1853delAA\nNM_000249.4:c.117-2del';
}

// Validate variant (single only, uses ferro directly)
async function validateVariant() {
    const isBatch = document.getElementById('mode-batch').classList.contains('active');

    if (isBatch) {
        alert('Validate is only available for single variants. Use Normalize for batch processing.');
        return;
    }

    const input = document.getElementById('hgvs-input').value.trim();

    if (!input) {
        alert('Please enter an HGVS variant');
        return;
    }

    showLoading();

    try {
        const response = await fetch('/api/v1/validate', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ hgvs: input })
        });

        const data = await response.json();

        if (response.ok) {
            displayValidateResult(data);
        } else {
            displayError(data);
        }
    } catch (error) {
        displayError({ error: 'request_failed', message: error.message });
    }
}

// Parse variant(s) - kept for backwards compatibility
async function parseVariant() {
    const isBatch = document.getElementById('mode-batch').classList.contains('active');

    if (isBatch) {
        await processBatch('parse');
    } else {
        await processSingle('parse');
    }
}

// Normalize variant(s)
async function normalizeVariant() {
    const isBatch = document.getElementById('mode-batch').classList.contains('active');

    if (isBatch) {
        await processBatch('normalize');
    } else {
        await processSingle('normalize');
    }
}

// Process single variant
async function processSingle(operation) {
    const input = document.getElementById('hgvs-input').value.trim();

    if (!input) {
        alert('Please enter an HGVS variant');
        return;
    }

    const tools = getSelectedTools();
    if (tools.length === 0) {
        alert('Please select at least one tool');
        return;
    }

    const errorMode = getErrorMode();

    showLoading();

    try {
        const response = await fetch(`/api/v1/${operation}`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                hgvs: input,
                tools: tools,
                error_mode: errorMode
            })
        });

        const data = await response.json();

        if (response.ok) {
            displaySingleResult(data);
        } else {
            displayError(data);
        }
    } catch (error) {
        displayError({ error: 'request_failed', message: error.message });
    }
}

// Process batch variants
async function processBatch(operation) {
    const input = document.getElementById('batch-input').value.trim();

    if (!input) {
        alert('Please enter HGVS variants');
        return;
    }

    const variants = input.split('\n').map(v => v.trim()).filter(v => v);

    if (variants.length === 0) {
        alert('Please enter at least one HGVS variant');
        return;
    }

    const tools = getSelectedTools();
    if (tools.length === 0) {
        alert('Please select at least one tool');
        return;
    }

    const errorMode = getErrorMode();

    showLoading();

    try {
        const response = await fetch(`/api/v1/batch/${operation}`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                variants: variants,
                tools: tools,
                error_mode: errorMode
            })
        });

        const data = await response.json();

        if (response.ok) {
            displayBatchResult(data);
        } else {
            displayError(data);
        }
    } catch (error) {
        displayError({ error: 'request_failed', message: error.message });
    }
}

// Show loading state
function showLoading() {
    const section = document.getElementById('results-section');
    const container = document.getElementById('results-container');

    section.style.display = 'block';
    container.innerHTML = '<div class="loading"><span class="spinner"></span> Processing...</div>';
    document.getElementById('processing-time').textContent = '';
}

// Display single result
function displaySingleResult(data) {
    const section = document.getElementById('results-section');
    const container = document.getElementById('results-container');
    const timeSpan = document.getElementById('processing-time');

    section.style.display = 'block';
    timeSpan.textContent = `${data.processing_time_ms}ms`;

    container.innerHTML = renderResultItem(data);
}

// Display batch result
function displayBatchResult(data) {
    const section = document.getElementById('results-section');
    const container = document.getElementById('results-container');
    const timeSpan = document.getElementById('processing-time');

    section.style.display = 'block';
    timeSpan.textContent = `${data.total_processing_time_ms}ms total`;

    let html = `<div class="batch-summary">
        <strong>Processed ${data.total_variants} variants</strong>
        (${data.successful_variants} successful)
    </div>`;

    for (const result of data.results) {
        html += renderResultItem(result);
    }

    container.innerHTML = html;
}

// Render a single result item
function renderResultItem(data) {
    const agreementClass = data.agreement.all_agree ? 'agree' : 'disagree';
    const agreementText = data.agreement.all_agree ? 'All Agree' : 'Disagree';

    // Determine success status class for color coding
    const totalTools = data.results.length;
    const successfulTools = data.agreement.successful_tools;
    let successClass = 'success-all';  // green
    if (successfulTools === 0) {
        successClass = 'success-none';  // red
    } else if (successfulTools < totalTools) {
        successClass = 'success-partial';  // yellow
    }

    let toolResultsHtml = '';
    let detailsHtml = '';

    for (const result of data.results) {
        const statusClass = result.success ? 'success' : 'error';

        // Check if output differs from input (normalized)
        const wasNormalized = result.success && result.output && result.output !== data.input;
        const normalizedIndicator = wasNormalized
            ? '<span class="normalized-badge">normalized</span>'
            : (result.success ? '<span class="unchanged-badge">unchanged</span>' : '');

        const outputHtml = result.success
            ? `<span class="tool-output">${escapeHtml(result.output)}</span>${normalizedIndicator}`
            : `<span class="tool-error">${escapeHtml(result.error)}</span>`;

        toolResultsHtml += `<div class="tool-result ${statusClass}">
            <span class="tool-name">${escapeHtml(result.tool)}</span>
            ${outputHtml}
            <span class="tool-time">${result.elapsed_ms}ms</span>
        </div>`;

        // Capture details from ferro tool if available
        if (result.tool === 'ferro' && result.details) {
            detailsHtml = renderComponentsTable(result.details);
        }
    }

    return `<div class="result-item">
        <div class="result-header">
            <span class="result-input">${escapeHtml(data.input)}</span>
            <div class="result-agreement">
                <span class="agreement-badge ${agreementClass}">${agreementText}</span>
                <span class="success-count ${successClass}">${data.agreement.successful_tools}/${data.results.length} tools succeeded</span>
            </div>
        </div>
        <div class="result-body">
            <div class="tool-results">
                ${toolResultsHtml}
            </div>
            ${detailsHtml}
        </div>
    </div>`;
}

// Display validate result
function displayValidateResult(data) {
    const section = document.getElementById('results-section');
    const container = document.getElementById('results-container');
    const timeSpan = document.getElementById('processing-time');

    section.style.display = 'block';
    timeSpan.textContent = `${data.processing_time_ms}ms`;

    const validClass = data.valid ? 'success' : 'error';
    const validText = data.valid ? 'Valid' : 'Invalid';
    const validBadgeClass = data.valid ? 'valid-badge' : 'invalid-badge';

    let html = `<div class="result-item">
        <div class="result-header">
            <span class="result-input">${escapeHtml(data.input)}</span>
            <span class="${validBadgeClass}">${validText}</span>
        </div>
        <div class="result-body">`;

    if (data.errors && data.errors.length > 0) {
        html += `<div class="validation-errors">
            <strong>Errors:</strong>
            <ul>`;
        for (const err of data.errors) {
            html += `<li>${escapeHtml(err)}</li>`;
        }
        html += `</ul></div>`;
    }

    if (data.components) {
        html += renderComponentsTable(data.components);
    }

    html += `</div></div>`;
    container.innerHTML = html;
}

// Render components/details table
function renderComponentsTable(components) {
    let html = `<div class="components-breakdown">
        <strong>Component Breakdown:</strong>
        <table class="components-table">
            <tr><td>Reference</td><td><code>${escapeHtml(components.reference)}</code></td></tr>
            <tr><td>Coordinate System</td><td><code>${escapeHtml(components.coordinate_system)}</code></td></tr>
            <tr><td>Variant Type</td><td><code>${escapeHtml(components.variant_type)}</code></td></tr>
            <tr><td>Position</td><td><code>${escapeHtml(components.position.display)}</code></td></tr>`;

    if (components.deleted) {
        html += `<tr><td>Deleted</td><td><code>${escapeHtml(components.deleted)}</code></td></tr>`;
    }
    if (components.inserted) {
        html += `<tr><td>Inserted</td><td><code>${escapeHtml(components.inserted)}</code></td></tr>`;
    }
    if (components.was_shifted !== undefined && components.was_shifted !== null) {
        html += `<tr><td>Was Shifted</td><td>${components.was_shifted ? 'Yes' : 'No'}</td></tr>`;
    }
    if (components.original_position) {
        html += `<tr><td>Original Position</td><td><code>${escapeHtml(components.original_position)}</code></td></tr>`;
    }

    html += `</table></div>`;
    return html;
}

// Display error with better formatting
function displayError(data) {
    const section = document.getElementById('results-section');
    const container = document.getElementById('results-container');

    section.style.display = 'block';
    document.getElementById('processing-time').textContent = '';

    // Map error codes to user-friendly messages
    const errorTitles = {
        'invalid_hgvs': 'Invalid HGVS Format',
        'bad_request': 'Invalid Request',
        'tool_unavailable': 'Tool Unavailable',
        'timeout': 'Request Timeout',
        'internal_error': 'Internal Error',
        'request_failed': 'Connection Error'
    };

    const errorTitle = errorTitles[data.error] || 'Error';
    const errorMessage = data.message || 'An unexpected error occurred.';

    // Add helpful hints based on error type
    let hint = '';
    if (data.error === 'invalid_hgvs') {
        hint = '<p class="error-hint">Check your HGVS syntax. Example: <code>NM_000249.4:c.350C>T</code></p>';
    } else if (data.error === 'request_failed') {
        hint = '<p class="error-hint">Unable to reach the server. Please check your connection and try again.</p>';
    }

    container.innerHTML = `<div class="result-item">
        <div class="result-header" style="background: var(--error);">
            <span style="color: white;">${escapeHtml(errorTitle)}</span>
        </div>
        <div class="result-body">
            <p>${escapeHtml(errorMessage)}</p>
            ${hint}
        </div>
    </div>`;
}

// Clear results
function clearResults() {
    document.getElementById('results-section').style.display = 'none';
    document.getElementById('results-container').innerHTML = '';
    document.getElementById('hgvs-input').value = '';
    document.getElementById('batch-input').value = '';
}

// Escape HTML to prevent XSS
function escapeHtml(text) {
    if (!text) return '';
    const div = document.createElement('div');
    div.textContent = text;
    return div.innerHTML;
}

// ==================== Operation Tab Navigation ====================

// Show operation panel
function showOperation(operation) {
    // Update tab buttons
    document.querySelectorAll('.op-tab-button').forEach(btn => {
        btn.classList.toggle('active', btn.dataset.op === operation);
    });

    // Update panels
    document.querySelectorAll('.operation-panel').forEach(panel => {
        panel.classList.toggle('active', panel.id === `op-${operation}`);
    });

    // Clear results when switching operations
    clearResults();
}

// ==================== Convert Operation ====================

// Set convert example
function setConvertExample(variant) {
    document.getElementById('convert-input').value = variant;
}

// Convert variant between coordinate systems
async function convertVariant() {
    const input = document.getElementById('convert-input').value.trim();

    if (!input) {
        alert('Please enter an HGVS variant');
        return;
    }

    const targetSystem = document.querySelector('input[name="target_system"]:checked').value;
    const includeAll = document.getElementById('include-all-conversions').checked;

    showLoading();

    try {
        const response = await fetch('/api/v1/convert', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                hgvs: input,
                target_system: targetSystem,
                include_all: includeAll
            })
        });

        const data = await response.json();

        if (response.ok) {
            displayConvertResult(data);
        } else {
            displayError(data);
        }
    } catch (error) {
        displayError({ error: 'request_failed', message: error.message });
    }
}

// Display convert result
function displayConvertResult(data) {
    const section = document.getElementById('results-section');
    const container = document.getElementById('results-container');
    const timeSpan = document.getElementById('processing-time');

    section.style.display = 'block';
    timeSpan.textContent = `${data.processing_time_ms}ms`;

    let html = `<div class="result-item">
        <div class="result-header">
            <span class="result-input">${escapeHtml(data.input)}</span>
            <span class="coord-badge">${escapeHtml(data.source_system)} → ${escapeHtml(data.target_system)}</span>
        </div>
        <div class="result-body">`;

    if (data.converted) {
        html += `<div class="conversion-result">
            <strong>Converted:</strong>
            <code>${escapeHtml(data.converted)}</code>
        </div>`;
    }

    if (data.all_conversions && data.all_conversions.length > 0) {
        html += `<div class="all-conversions">
            <strong>All Conversions:</strong>
            <ul>`;
        for (const conv of data.all_conversions) {
            html += `<li><code>${escapeHtml(conv.hgvs)}</code> (${escapeHtml(conv.system)})</li>`;
        }
        html += `</ul></div>`;
    }

    if (data.error) {
        html += `<div class="conversion-error">
            <strong>Note:</strong> ${escapeHtml(data.error)}
        </div>`;
    }

    html += `</div></div>`;
    container.innerHTML = html;
}

// ==================== Effect Operation ====================

// Set effect example
function setEffectExample(variant) {
    document.getElementById('effect-input').value = variant;
}

// Predict effect
async function predictEffect() {
    const input = document.getElementById('effect-input').value.trim();

    if (!input) {
        alert('Please enter an HGVS variant');
        return;
    }

    const includeNmd = document.getElementById('include-nmd').checked;

    showLoading();

    try {
        const response = await fetch('/api/v1/effect', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                hgvs: input,
                include_nmd: includeNmd
            })
        });

        const data = await response.json();

        if (response.ok) {
            displayEffectResult(data);
        } else {
            displayError(data);
        }
    } catch (error) {
        displayError({ error: 'request_failed', message: error.message });
    }
}

// Display effect result
function displayEffectResult(data) {
    const section = document.getElementById('results-section');
    const container = document.getElementById('results-container');
    const timeSpan = document.getElementById('processing-time');

    section.style.display = 'block';
    timeSpan.textContent = `${data.processing_time_ms}ms`;

    let html = `<div class="result-item">
        <div class="result-header">
            <span class="result-input">${escapeHtml(data.input)}</span>
        </div>
        <div class="result-body">`;

    if (data.effect) {
        const impactClass = data.effect.impact.toLowerCase();
        html += `<div class="effect-result">
            <table class="components-table">
                <tr><td>SO Term</td><td><code>${escapeHtml(data.effect.so_term)}</code></td></tr>
                <tr><td>Effect</td><td><strong>${escapeHtml(data.effect.name)}</strong></td></tr>
                <tr><td>Description</td><td>${escapeHtml(data.effect.description)}</td></tr>
                <tr><td>Impact</td><td><span class="impact-badge impact-${impactClass}">${escapeHtml(data.effect.impact)}</span></td></tr>
            </table>
        </div>`;
    }

    if (data.protein_consequence) {
        html += `<div class="protein-result">
            <strong>Protein Consequence:</strong>
            <code>${escapeHtml(data.protein_consequence.hgvs_p)}</code>
        </div>`;
    }

    if (data.nmd_prediction) {
        html += `<div class="nmd-result">
            <strong>NMD Prediction:</strong>
            ${data.nmd_prediction.predicted ? 'Yes' : 'No'}
            (confidence: ${(data.nmd_prediction.confidence * 100).toFixed(0)}%)
            <br><small>${escapeHtml(data.nmd_prediction.reason)}</small>
        </div>`;
    }

    if (data.error) {
        html += `<div class="conversion-error">
            <strong>Note:</strong> ${escapeHtml(data.error)}
        </div>`;
    }

    html += `</div></div>`;
    container.innerHTML = html;
}

// ==================== Liftover Operation ====================

// Set liftover example
function setLiftoverExample(position) {
    document.getElementById('liftover-input').value = position;
}

// Liftover position
async function liftoverPosition() {
    const input = document.getElementById('liftover-input').value.trim();

    if (!input) {
        alert('Please enter a genomic position');
        return;
    }

    const fromBuild = document.getElementById('from-build').value;
    const toBuild = document.getElementById('to-build').value;

    if (fromBuild === toBuild) {
        alert('Source and target builds must be different');
        return;
    }

    showLoading();

    try {
        const response = await fetch('/api/v1/liftover', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                position: input,
                from_build: fromBuild,
                to_build: toBuild
            })
        });

        const data = await response.json();

        if (response.ok) {
            displayLiftoverResult(data);
        } else {
            displayError(data);
        }
    } catch (error) {
        displayError({ error: 'request_failed', message: error.message });
    }
}

// Display liftover result
function displayLiftoverResult(data) {
    const section = document.getElementById('results-section');
    const container = document.getElementById('results-container');
    const timeSpan = document.getElementById('processing-time');

    section.style.display = 'block';
    timeSpan.textContent = `${data.processing_time_ms}ms`;

    let html = `<div class="result-item">
        <div class="result-header">
            <span class="result-input">${escapeHtml(data.input)}</span>
            <span class="coord-badge">${escapeHtml(data.from_build)} → ${escapeHtml(data.to_build)}</span>
        </div>
        <div class="result-body">`;

    if (data.converted) {
        html += `<div class="conversion-result">
            <strong>Converted Position:</strong>
            <code>${escapeHtml(data.converted)}</code>
        </div>`;
    }

    if (data.hgvs_g) {
        html += `<div class="conversion-result">
            <strong>HGVS Genomic:</strong>
            <code>${escapeHtml(data.hgvs_g)}</code>
        </div>`;
    }

    if (data.chain_region) {
        html += `<div class="chain-region">
            <strong>Region:</strong>
            <code>${escapeHtml(data.chain_region)}</code>
        </div>`;
    }

    if (data.error) {
        html += `<div class="conversion-error">
            <strong>Note:</strong> ${escapeHtml(data.error)}
        </div>`;
    }

    html += `</div></div>`;
    container.innerHTML = html;
}

// ==================== VCF Operation ====================

// Show VCF mode
function showVcfMode(mode) {
    document.querySelectorAll('.vcf-mode-tabs .tab-button').forEach(btn => {
        btn.classList.toggle('active', btn.textContent.toLowerCase().includes(mode.replace('-', ' ')));
    });

    document.querySelectorAll('.vcf-mode').forEach(panel => {
        panel.classList.toggle('active', panel.id === `vcf-mode-${mode}`);
    });
}

// Set HGVS to VCF example
function setHgvsToVcfExample(hgvs) {
    document.getElementById('hgvs-to-vcf-input').value = hgvs;
}

// Set VCF to HGVS example
function setVcfToHgvsExample(chrom, pos, ref, alt) {
    document.getElementById('vcf-chrom').value = chrom;
    document.getElementById('vcf-pos').value = pos;
    document.getElementById('vcf-ref').value = ref;
    document.getElementById('vcf-alt').value = alt;
}

// HGVS to VCF conversion
async function hgvsToVcf() {
    const input = document.getElementById('hgvs-to-vcf-input').value.trim();

    if (!input) {
        alert('Please enter an HGVS variant');
        return;
    }

    const build = document.getElementById('hgvs-to-vcf-build').value;

    showLoading();

    try {
        const response = await fetch('/api/v1/hgvs-to-vcf', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                hgvs: input,
                build: build
            })
        });

        const data = await response.json();

        if (response.ok) {
            displayHgvsToVcfResult(data);
        } else {
            displayError(data);
        }
    } catch (error) {
        displayError({ error: 'request_failed', message: error.message });
    }
}

// VCF to HGVS conversion
async function vcfToHgvs() {
    const chrom = document.getElementById('vcf-chrom').value.trim();
    const pos = document.getElementById('vcf-pos').value.trim();
    const ref = document.getElementById('vcf-ref').value.trim();
    const alt = document.getElementById('vcf-alt').value.trim();

    if (!chrom || !pos || !ref || !alt) {
        alert('Please fill in all VCF fields');
        return;
    }

    const build = document.getElementById('vcf-to-hgvs-build').value;
    const transcript = document.getElementById('vcf-transcript').value.trim() || null;

    showLoading();

    try {
        const response = await fetch('/api/v1/vcf-to-hgvs', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                chrom: chrom,
                pos: parseInt(pos, 10),
                ref: ref,
                alt: alt,
                build: build,
                transcript: transcript
            })
        });

        const data = await response.json();

        if (response.ok) {
            displayVcfToHgvsResult(data);
        } else {
            displayError(data);
        }
    } catch (error) {
        displayError({ error: 'request_failed', message: error.message });
    }
}

// Display HGVS to VCF result
function displayHgvsToVcfResult(data) {
    const section = document.getElementById('results-section');
    const container = document.getElementById('results-container');
    const timeSpan = document.getElementById('processing-time');

    section.style.display = 'block';
    timeSpan.textContent = `${data.processing_time_ms}ms`;

    let html = `<div class="result-item">
        <div class="result-header">
            <span class="result-input">${escapeHtml(data.input)}</span>
            <span class="coord-badge">HGVS → VCF</span>
        </div>
        <div class="result-body">`;

    if (data.vcf) {
        html += `<div class="vcf-result">
            <table class="components-table">
                <tr><td>CHROM</td><td><code>${escapeHtml(data.vcf.chrom)}</code></td></tr>
                <tr><td>POS</td><td><code>${data.vcf.pos}</code></td></tr>
                <tr><td>REF</td><td><code>${escapeHtml(data.vcf.ref)}</code></td></tr>
                <tr><td>ALT</td><td><code>${escapeHtml(data.vcf.alt)}</code></td></tr>
                <tr><td>Build</td><td>${escapeHtml(data.vcf.build)}</td></tr>
            </table>
        </div>`;
    }

    if (data.error) {
        html += `<div class="conversion-error">
            <strong>Note:</strong> ${escapeHtml(data.error)}
        </div>`;
    }

    html += `</div></div>`;
    container.innerHTML = html;
}

// Display VCF to HGVS result
function displayVcfToHgvsResult(data) {
    const section = document.getElementById('results-section');
    const container = document.getElementById('results-container');
    const timeSpan = document.getElementById('processing-time');

    section.style.display = 'block';
    timeSpan.textContent = `${data.processing_time_ms}ms`;

    const vcf = data.vcf;
    let html = `<div class="result-item">
        <div class="result-header">
            <span class="result-input">${escapeHtml(vcf.chrom)}:${vcf.pos} ${escapeHtml(vcf.ref)}>${escapeHtml(vcf.alt)}</span>
            <span class="coord-badge">VCF → HGVS</span>
        </div>
        <div class="result-body">`;

    if (data.hgvs_g) {
        html += `<div class="conversion-result">
            <strong>Genomic (g.):</strong>
            <code>${escapeHtml(data.hgvs_g)}</code>
        </div>`;
    }

    if (data.hgvs_c) {
        html += `<div class="conversion-result">
            <strong>Coding (c.):</strong>
            <code>${escapeHtml(data.hgvs_c)}</code>
        </div>`;
    }

    if (data.hgvs_p) {
        html += `<div class="conversion-result">
            <strong>Protein (p.):</strong>
            <code>${escapeHtml(data.hgvs_p)}</code>
        </div>`;
    }

    if (data.error) {
        html += `<div class="conversion-error">
            <strong>Note:</strong> ${escapeHtml(data.error)}
        </div>`;
    }

    html += `</div></div>`;
    container.innerHTML = html;
}
