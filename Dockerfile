FROM python:3.8.4
LABEL maintainer="Jaeyoung Chun (chunj@mskcc.org)"

RUN pip install scanpy==1.5.1

COPY ./scripts/archr2adata.py /opt/archr2adata.py