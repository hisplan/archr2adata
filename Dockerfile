FROM python:3.8.4
LABEL maintainer="Jaeyoung Chun (chunj@mskcc.org)"

RUN pip install --upgrade pip
RUN pip install scanpy==1.8.2

COPY ./scripts/archr2adata.py /opt/archr2adata.py
