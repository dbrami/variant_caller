FROM python:3.9-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY src/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy source code
COPY src/variant_caller /app/variant_caller
COPY src/setup.py .

# Install the package
RUN pip install -e .

# Set entrypoint to Python
ENTRYPOINT ["python", "-m"]
