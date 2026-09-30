FROM python:3.12-slim

LABEL org.opencontainers.image.title="OncoFlow Atlas" \
      org.opencontainers.image.description="Reproducible cancer genomics multi-omic workflow" \
      org.opencontainers.image.version="1.0.0"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /opt/oncoflow
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
COPY . .

CMD ["python", "bin/analyse_multiomics.py", "--help"]
