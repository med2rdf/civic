# setting base image
FROM continuumio/anaconda:latest

ADD tsv2rdf_civic.py /tmp
ADD templ_civic.ttl /tmp
ADD templ_civic.ttl.prefix /tmp
ADD templ_civic.ttl.variant /tmp
ADD templ_civic.ttl.evi /tmp
ADD templ_civic.ttl.drug /tmp
ADD templ_civic.ttl.disease /tmp
ADD so_table.tsv /tmp
ADD drug_table.tsv /tmp
ADD disease_table.tsv /tmp
ENTRYPOINT ["python", "/tmp/tsv2rdf_civic.py"]
CMD ["-c","/mnt/tsv2rdf_civic.json"]
