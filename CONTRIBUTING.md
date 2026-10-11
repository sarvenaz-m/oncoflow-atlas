# Contributing

Keep input provenance, mathematical definitions and rendered results consistent. Edit the browser sources in `src/oncoflow_atlas/web/` and rebuild `docs/index.html`; do not patch only the generated HTML.

For analysis changes, run the regression suite and independent numerical/parity checks. For interface changes, run the browser checks. Follow the commands in the README and update the recorded evidence when intentionally rebuilding release artifacts. An implementation change must not be presented as new biological evidence.

Keep the pinned upstream source files unchanged. To introduce a new study or source revision, create a separate lock/data card with specimen selection, feature exclusions, license and units before interpreting results. Never replace missing clinical information with invented values.

The public data terms are separate from the software license. Preserve source attribution and the complete open adapted dataset when redistributing data-bearing reports.
