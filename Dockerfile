# Dockerfile for hfe-variant-prevalence

FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Copy source code
COPY . .

# Install the adapter after copying the package and its README.
RUN python -m pip install --only-binary=:all: --no-cache-dir .

# Create logs directory
RUN mkdir -p logs

# Set environment variables
ENV PYTHONPATH=/app/src
ENV PYTHONUNBUFFERED=1

# Default command
CMD ["python", "create_knowledge_graph.py"]
