FROM python:3.10-slim-buster

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libboost-dev \
    libboost-system-dev \
    libboost-thread-dev \
    libeigen3-dev \
    libxrender1 \
    libsm6 \
    libxext6 \
    libgl1-mesa-glx \
    && rm -rf /var/lib/apt/lists/*

## Set the working directory
WORKDIR /app

# Copy dependency files
COPY requirements.txt requirements.txt /app/

# Install Poetry
RUN pip install -r requirements.txt


# Copy application code
COPY /app/Standardization_Pipeline /app/Standardization_Pipeline
COPY /app/function_app.py /app/function_app.py

# Expose port 8000
EXPOSE 8000

# Command to run the app
CMD ["uvicorn", "function_app:app", "--host", "0.0.0.0", "--port", "8000"]
