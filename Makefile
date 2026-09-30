.PHONY: demo nextflow test lint quality clean

demo:
	python bin/generate_demo.py --outdir data/demo --samples 36 --seed 19
	python bin/analyse_multiomics.py --expression data/demo/expression.csv --variants data/demo/variants.csv --metadata data/demo/metadata.csv --outdir results/demo --cohort-name demo-cohort --min-vaf 0.05 --top-n 12

nextflow:
	nextflow run main.nf -profile local -resume

test:
	pytest -q

lint:
	ruff check bin tests

quality: lint test

clean:
	rm -rf results/demo results/nextflow .pytest_cache .ruff_cache

showcase: demo
	cp results/demo/report.html docs/demo.html
	@echo "Showcase refreshed: docs/demo.html"
