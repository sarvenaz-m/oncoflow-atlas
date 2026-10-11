FROM python:3.12-slim
LABEL org.opencontainers.image.title="OncoFlow Atlas" \
      org.opencontainers.image.version="2.0.0"
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
WORKDIR /opt/oncoflow
COPY requirements-lock.txt ./
RUN pip install --no-cache-dir -r requirements-lock.txt
COPY pyproject.toml ./
COPY src ./src
RUN pip install --no-cache-dir --no-deps .
CMD ["oncoflow-atlas", "--help"]
