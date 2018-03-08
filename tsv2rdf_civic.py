#coding: UTF-8
'''
Created on 2018/03/04
'''
import sys,os
from jinja2 import Template, Environment, FileSystemLoader
from optparse import OptionParser
import json
import numpy as np
import logging.handlers
import logging
import pandas as pd
import csv
import re

RDF_TEMPLATE='templin_civic.ttl'
EXEC_PATH = '/tmp'
MOUNT_PATH = '/mnt'
TEMPLATE_BODY = 'templ_civic.ttl'
TEMPLATE_EVI = 'templ_civic.ttl.evi'
TEMPLATE_PREFIX = 'templ_civic.ttl.prefix'
TEMPLATE_VARIANT = 'templ_civic.ttl.variant'
TEMPLATE_DISEASE = 'templ_civic.ttl.disease'
TEMPLATE_DRUG = 'templ_civic.ttl.drug'

class Gene(object):
    def __init__(self, gene_id, gene, entrez_id, data ):
        '''
        @summary: Store gene data (1 row) in the Message object
        '''

        self.gene_id = gene_id
        self.gene = gene
        self.entrez_id = entrez_id
        
        #Delete duplicate variant_id
        data_uniq = data.drop_duplicates(subset=['gene_id', 'variant_id'], inplace=False )

        self.Variants = [ Variants(comma, variant_id)
                             for comma, variant_id in enumerate(data_uniq['variant_id'])] 
        
        return

class Variants(object):
    def __init__(self, comma, variant_id ):
        self.variant_id = variant_id
        self.comma = comma      # 0: No comma >1: comma present

        return

class Variant(object):
    def __init__(self, variant_id, variant_row, evidence_rows, so_table_dic, disease_table_dic, sequence_no=0 ):
        '''
        @param evidence_rows: Evidence DataFrame
        @summary: Store variant data (1 row) in the Message object
        '''
        self.has_reference_bases = False
        self.has_variant_bases = False
        self.has_hgvs_expressions = False
        self.has_reference_build = False
        self.has_representative_transcript = False
        self.has_ensembl_version = False
        self.has_variant_origin = False
        self.has_variant_summary = False
        self.has_start = False
        self.has_gene_id = False
        self.has_doid = False

        self.variant_id = variant_id
        self.variant = variant_row.variant.values[0].decode('utf-8')

        self.variant_types=variant_row.variant_types.values[0]
        self.hgvs_expressions = variant_row.hgvs_expressions.values[0]

        list_variant_types = map(lambda x: so_table_dic[x]  , self.variant_types.split(','))

        self.VariantTypes = [ VariantType(comma, variant_type)
                    for comma, variant_type in enumerate(list_variant_types)] 

        if evidence_rows.empty:
            self.has_evidence = False
            if len(self.hgvs_expressions) != 0:self.has_hgvs_expressions = True
            #logger.info("Variant has no evidence: variant_id=%s" % (self.variant_id))
        else:
            self.has_evidence = True

            self.Publications = [ PublicationId(comma, row) for comma, row in evidence_rows.iterrows()]

            #Delete duplicate rows
            evidence_row = evidence_rows.drop_duplicates(subset=['variant_id'], inplace=False )
            evidence_row = evidence_row.reset_index(drop=True)    #Reset index(0 origin)

            self.disease_code = disease_table_dic[evidence_row.disease.values[0]]
            self.variant_origin=evidence_row.variant_origin.values[0]      
            self.variant_summary=evidence_row.variant_summary.values[0].replace('"','\\"')
            
            self.start=formatter(evidence_row.start.values[0])
            self.stop=formatter(evidence_row.stop.values[0])

 
            if len(str(evidence_row.start2.values[0])) != 0:
                self.has_start2 = True
                self.start2=formatter(evidence_row.start2.values[0])
                self.stop2=formatter(evidence_row.stop2.values[0])
                self.chromosome2=evidence_row.chromosome2.values[0]
                self.sequence_no_1 = str(sequence_no) + '_1'
                self.sequence_no_2 = str(sequence_no) + '_2'
            else:
                self.has_start2 = False
            
            self.chromosome=evidence_row.chromosome.values[0]
            self.gene_id=evidence_row.gene_id.values[0]

            self.variant_origin=evidence_row.variant_origin.values[0]
            
            variant_summary = evidence_row.variant_summary.values[0].replace('"','\\"')
            
            if isinstance(variant_summary, str):
                self.variant_summary=variant_summary.decode('utf-8')
            else:
                self.variant_summary=''

            self.reference_bases=evidence_row.reference_bases.values[0]
            self.variant_bases=evidence_row.variant_bases.values[0]
            self.reference_build=evidence_row.reference_build.values[0]
            self.representative_transcript=evidence_row.representative_transcript.values[0]
            self.ensembl_version=formatter(evidence_row.ensembl_version.values[0])
            wk_df = evidence_rows.drop_duplicates(subset=['doid'], inplace=False )
            #wk_df = wk_df[~wk_df['doid'].isnull()].reset_index(drop=True)
            
            self.Doids = []
            comma = 0
            for doid in evidence_rows.drop_duplicates(subset=['doid'], inplace=False )['doid']:
                if len(str(doid)) != 0: 
                    self.Doids.append(Doid(comma, doid))
                    comma +=1
                else:
                    pass

            if len(self.Doids) != 0:self.has_doid = True
            if len(self.variant_types) != 0:self.has_variant_types = True
            if len(self.reference_bases) != 0:self.has_reference_bases = True
            if len(self.variant_bases) != 0:self.has_variant_bases = True
            if len(self.reference_build) != 0:self.has_reference_build = True
            if len(self.representative_transcript) != 0:self.has_representative_transcript = True
            if len(self.ensembl_version) != 0:self.has_ensembl_version = True
            if len(self.variant_origin) != 0:self.has_variant_origin = True
            if len(self.variant_summary) != 0:self.has_variant_summary = True
            if len(self.start) != 0:self.has_start = True
            if len(str(self.gene_id)) != 0:self.has_gene_id = True

        return

class VariantType(object):
    def __init__(self, comma, variant_type):
        '''
        @param comma: Variables for controlling output of commas when outputting multiple publications
        @summary: Store multiple VariantType in the VariantTypes object
        '''
        self.comma = comma      # 0: No comma >1: comma present
        self.variant_type = variant_type

        return

class PublicationId(object):
    def __init__(self, comma, row):
        '''
        @param comma: Variables for controlling output of commas when outputting multiple publications
        @summary: Store multiple evidence_id and doid in the Publications object
        '''
        self.has_doid = False

        self.comma = comma      # 0: No comma >1: Comma present
        self.evidence_id=row.evidence_id
        self.doid=formatter(row.doid)
        
        if len(self.doid) != 0:
            self.has_doid = True

        return

class Publication(object):
    def __init__(self, row, drug_table_dic, drug_label_dic):

        self.has_rating = False
        self.has_evidence_direction = False
        self.has_clinical_significance = False
        self.has_drugs = False
        self.has_doid = False

        self.evidence_id=row.evidence_id
        self.pubmed_id=row.pubmed_id
        self.citation=row.citation.decode('utf-8')
        self.rating=row.rating
        self.doid=formatter(row.doid)

        self.evidence_status=row.evidence_status
        self.evidence_type=row.evidence_type
        self.evidence_direction=row.evidence_direction
        self.clinical_significance=row.clinical_significance
        self.evidence_level=row.evidence_level
        self.drug_description=row.drugs


        #Split drugs
        wk_drugs_uniq = get_drug_list(self.drug_description, drug_table_dic)
        self.Drugs =  [ Drug(comma, drugname, drug_label_dic[drugname]) for comma, drugname in enumerate(wk_drugs_uniq)]

        self.evidence_statement=row.evidence_statement.replace('"','\\"').decode('utf-8')
        self.disease=row.disease.decode('utf-8')

        if len(self.doid) != 0:self.has_doid = True
        if len(self.drug_description)>0:self.has_drugs = True
        if len(str(self.rating)) != 0:self.has_rating = True
        if len(self.evidence_direction) != 0:self.has_evidence_direction = True
        if len(self.clinical_significance) != 0:self.has_clinical_significance = True

        return

class Doid(object):
    def __init__(self, comma, doid):
        self.has_doid = False

        self.comma = comma      # 0: No comma >1: Comma present
        self.doid=formatter(doid)

        if len(str(self.doid)) != 0:self.has_doid = True

        return

class Disease(object):
    def __init__(self, row, disease_table_dic):
        self.has_doid = False

        self.disease_code = disease_table_dic[row.disease]
        self.doid = formatter(row.doid)
        self.disease = row.disease.decode('utf-8')

        if len(self.doid) != 0:self.has_doid = True

        return

    
class Drug(object):
    def __init__(self, comma, name, label):

        self.comma = comma      # 0: No comma >1: Comma present
        self.name=name  
        self.label=label

        return

def formatter(val):
    retval=None
    
    if isinstance(val, float):
        retval = '{:.0f}'.format(val)
    else:
        retval = ''

    return retval

def get_drug_list(drugs, drug_table_dic, name_type='StandardName'):
    '''
    @param drugs: original drugs
    @param drug_table_dic: drug dictionary
    '''

    wk_drugs = []

    for k,v in drug_table_dic.items():
        regex = '.*%s.*' % k.strip().replace('(', '\(').replace(')', '\)').replace('[', '\[').replace(']', '\]')
        m = re.match(regex, drugs)
        if m:
            wk_drugs.append(v[name_type])

    drugs_uniq = list(set(wk_drugs))

    return drugs_uniq

def main(fw, template, evidence_summaries, variant_summaries, conversion_tables ):

    #SequenceOntology table
    so_table_dic = {}
    reader = csv.DictReader(open(conversion_tables['sequence_ontology'],'r'), delimiter="\t")
    for row in reader:
        so_table_dic[row['Civic']] = row['SequenceOntology']

    # Disease table
    disease_table_dic = {}
    reader = csv.DictReader(open(conversion_tables['disease'],'r'), delimiter="\t")
    for row in reader:
        disease_table_dic[row['disease']] = row['URI']

    variantDF=pd.read_csv(variant_summaries, delimiter='\t')
    variantDF.fillna('', inplace=True)
    evidenceDF=pd.read_csv(evidence_summaries, delimiter='\t')
    evidenceDF.set_index('variant_id')
    evidenceDF.fillna('', inplace=True)

    '''
    Processing gene-template (gene-body)
    '''
    env = Environment(loader = FileSystemLoader('.', encoding='utf8'), autoescape = False)
    imageTemplate = env.get_template(template['body'])

    geneDf = evidenceDF.drop_duplicates(['gene_id','gene'])
        
    for i, row in geneDf.iterrows():

        evidenceDF_selected = evidenceDF[ evidenceDF["gene_id"] == row["gene_id"]]

        #Message Body
        msgObject = Gene(row["gene_id"], row["gene"], row["entrez_id"], evidenceDF_selected)

        namespace = dict(message=msgObject)
        FeedContent = imageTemplate.render(namespace).encode('utf-8')

        fw.write(FeedContent)

    '''
    Processing variant-template
    '''
    imageTemplate = env.get_template(template['variant'])

    for i, row in variantDF.iterrows():

        variant_row = variantDF[ variantDF["variant_id"] == row["variant_id"]]
        evidence_rows = evidenceDF[ evidenceDF["variant_id"] == row["variant_id"]]
        evidence_rows = evidence_rows.reset_index(drop=True)    #Reset index(0 origin)

        msgObject = Variant(row["variant_id"], variant_row, evidence_rows, so_table_dic, disease_table_dic, sequence_no=i)

        namespace = dict(message=msgObject)
        FeedContent = imageTemplate.render(namespace).encode('utf-8')

        fw.write(FeedContent)

    '''
    Processing evidence-template
    '''
    #Drug table
    drug_table_dic = {}
    drug_label_dic = {}
    reader = csv.DictReader(open(conversion_tables['drugs'],'r'), delimiter="\t")
    for row in reader:
        drug_table_dic[row['DrugName']] = dict( (('StandardName',row['StandardName']), ('DrugLabel',row['DrugLabel'])) )

    for k, v in drug_table_dic.items():
        drug_label_dic[v['StandardName']] = v['DrugLabel']

    imageTemplate = env.get_template(template['evidence'])

    #evidenceDF_selected = evidenceDF.drop_duplicates(['doid'])
    evidenceDF_selected = evidenceDF

    for i, row in evidenceDF_selected.iterrows():

        msgObject = Publication(row, drug_table_dic, drug_label_dic)

        namespace = dict(message=msgObject)
        FeedContent = imageTemplate.render(namespace).encode('utf-8')
        fw.write(FeedContent)

    '''
    Processing disease-template
    '''       
    imageTemplate = env.get_template(template['disease'])

    evidenceDF_selected = evidenceDF.drop_duplicates(['disease'])
    for i, row in evidenceDF_selected.iterrows():
        msgObject = Disease(row, disease_table_dic)

        namespace = dict(message=msgObject)
        FeedContent = imageTemplate.render(namespace).encode('utf-8')
        fw.write(FeedContent)
    
    '''
    Processing drug-template
    '''
    imageTemplate = env.get_template(template['drug'])

    drugs_list = []
    evidenceDF_selected = evidenceDF.drop_duplicates(['drugs'])
    for i, row in evidenceDF_selected.iterrows():
        drugs_list.extend(get_drug_list(row.drugs, drug_table_dic, name_type='StandardName'))        

    drugs_list_uniq = list(set(drugs_list))

    msgObject = [ Drug(0, drug, drug_label_dic[drug]) for drug in drugs_list_uniq ]
    namespace = dict(message=msgObject)
    FeedContent = imageTemplate.render(namespace).encode('utf-8')
    fw.write(FeedContent)
        
    return

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option( "-c", type="string", help="config_file", dest="config", default= 'tsv2rdf_civic.json')
    (options, args) = parser.parse_args()

    #Log setting
    argvs = sys.argv
    LOG_FILENAME = '%s.log' % argvs[0].split('.')[0]
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s',level=logging.DEBUG)
    logger = logging.getLogger()

    #File log
    file_log = logging.handlers.RotatingFileHandler(filename=LOG_FILENAME)
    file_log.setLevel(logging.INFO)
    file_log.setFormatter(logging.Formatter('%(asctime)s;%(levelname)s:%(message)s'))

    #Log handler setting
    logging.getLogger().addHandler(file_log)

    #Reading the configuration file
    try:
        f = open(options.config ,"r")
    except Exception as e:
        logger.error("%s: %s" % (e.message,options.config))
        sys.exit()
    
    f = open(options.config ,"r")
    config = json.load(f)

    output_file = config['output_file']
    template = config['template']
    conversion_tables = config['conversion_tables']
    data_path = config['data_path']

    #Output file
    if os.getcwd() == '/tmp':
        output_file = os.path.join(MOUNT_PATH, output_file)
        fw = open(output_file, 'w')
    else:
        fw = open(output_file, 'w')

    #Check template file
    if os.path.exists(template['body']) != True:
        template['body'] = os.path.join(EXEC_PATH, template['body'])
        logger.info("Set template file path: %s" % template['body'])
    if os.path.exists(template['evidence']) != True:
        template['evidence'] = os.path.join(EXEC_PATH, template['evidence'])
        logger.info("Set template file path: %s" % template['evidence'])
    if os.path.exists(template['prefix']) != True:
        template['prefix'] = os.path.join(EXEC_PATH, template['prefix'])
        logger.info("Set template file path: %s" % template['prefix'])
    if os.path.exists(template['variant']) != True:
        template['variant'] = os.path.join(EXEC_PATH, template['variant'])
        logger.info("Set template file path: %s" % template['variant'])
    if os.path.exists(template['disease']) != True:
        template['disease'] = os.path.join(EXEC_PATH, template['disease'])
        logger.info("Set template file path: %s" % template['disease'])
    if os.path.exists(template['drug']) != True:
        template['drug'] = os.path.join(EXEC_PATH, template['drug'])
        logger.info("Set template file path: %s" % template['drug'])
        
    if os.path.exists(conversion_tables['sequence_ontology']) != True:
        conversion_tables['sequence_ontology'] = os.path.join(EXEC_PATH, conversion_tables['sequence_ontology'])
        logger.info("Set conversion_tables file path: %s" % conversion_tables['sequence_ontology'])
    if os.path.exists(conversion_tables['drugs']) != True:
        conversion_tables['drugs'] = os.path.join(EXEC_PATH, conversion_tables['drugs'])
        logger.info("Set conversion_tables file path: %s" % conversion_tables['drugs'])
    if os.path.exists(conversion_tables['disease']) != True:
        conversion_tables['disease'] = os.path.join(EXEC_PATH, conversion_tables['disease'])
        logger.info("Set conversion_tables file path: %s" % conversion_tables['disease'])

    if os.path.exists(data_path) != True:
        data_path = os.path.join(MOUNT_PATH, data_path)
        logger.info("Set data_path: %s" % data_path)

    #Output prefix
    fp = open(template['prefix'], 'r')
    prefixes = fp.readlines()
    for prefix in prefixes:
        fw.write(prefix)

    evidence_summaries=os.path.join(data_path, config['input_file']['evidence_summaries'])
    variant_summaries=os.path.join(data_path, config['input_file']['variant_summaries'])

    main( fw, template, evidence_summaries, variant_summaries, conversion_tables )
    pass
