FROM python:3.9-slim

WORKDIR /app

# Install system dependencies including procps and tabix
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    procps \
    tabix \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY src/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy source code
COPY src/variant_caller /app/variant_caller
COPY src/setup.py .

# Install the package in development mode
RUN pip install -e . && \
    pip list && \
    which detect_sv && \
    detect_sv --help

# Set working directory for input/output
WORKDIR /data