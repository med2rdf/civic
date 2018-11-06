# setting base image
FROM continuumio/anaconda3:latest

ADD tsv2rdf_civic.py /tmp
ADD check_tsv.py /tmp
ADD templ_civic.ttl /tmp
ADD templ_civic.ttl.prefix /tmp
ADD templ_civic.ttl.variant /tmp
ADD templ_civic.ttl.evi /tmp
ADD templ_civic.ttl.drug /tmp
ADD templ_civic.ttl.disease /tmp
ENTRYPOINT ["python3", "/tmp/tsv2rdf_civic.py"]
CMD ["-c","/mnt/tsv2rdf_civic.json"]
