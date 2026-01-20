#!/usr/bin/env python3
"""
NCBI Entrez CLI tool for querying clinical metagenomics samples
with both short and long read data available.

Enhanced with BioProject search, PubMed links, and validation mode.
"""

import argparse
import sys
import yaml
import logging
from typing import List, Optional, Dict, Tuple
from urllib.parse import quote
from urllib.request import urlopen
from urllib.error import HTTPError, URLError
import xml.etree.ElementTree as ET
import time
import json
import re

# Configure logging
logger = logging.getLogger('entrez_tool')
logger.setLevel(logging.INFO)

class EntrezQueryTool:
    """Tool for querying NCBI Entrez databases for metagenomics samples."""
    
    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    def __init__(self, email: str = "user@example.com", api_key: Optional[str] = None):
        self.email = email
        self.api_key = api_key
        self.delay = 0.34 if not api_key else 0.1
    
    def _build_url(self, endpoint: str, params: Dict[str, str]) -> str:
        """Build Entrez API URL with parameters."""
        params['email'] = self.email
        if self.api_key:
            params['api_key'] = self.api_key
        
        param_str = '&'.join([f"{k}={quote(str(v))}" for k, v in params.items()])
        return f"{self.BASE_URL}{endpoint}?{param_str}"
    
    def _make_request(self, url: str) -> str:
        """Make HTTP request with rate limiting."""
        time.sleep(self.delay)
        try:
            with urlopen(url, timeout=30) as response:
                return response.read().decode('utf-8')
        except HTTPError as e:
            logger.error(f"HTTP Error {e.code}: {e.reason}")
            return None
        except URLError as e:
            logger.error(f"URL Error: {e.reason}")
            return None
    
    def search_pubmed(self, query: str, retmax: int = 20) -> List[Dict]:
        """Search PubMed for publications and return article details with links."""
        params = {
            'db': 'pubmed',
            'term': query,
            'retmax': str(retmax),
            'retmode': 'json'
        }
        
        url = self._build_url('esearch.fcgi', params)
        logger.info(f"[PubMed Search] Query: {query}")
        
        response = self._make_request(url)
        if not response:
            return []
        
        try:
            data = json.loads(response)
            id_list = data.get('esearchresult', {}).get('idlist', [])
            count = data.get('esearchresult', {}).get('count', '0')
            logger.info(f"[PubMed] Found {count} publications, retrieving {len(id_list)}")
            
            if not id_list:
                return []
            
            # Fetch publication details
            return self._fetch_pubmed_details(id_list)
        
        except json.JSONDecodeError:
            logger.error("Error parsing PubMed results")
            return []
    
    def _fetch_pubmed_details(self, pmid_list: List[str]) -> List[Dict]:
        """Fetch PubMed article details."""
        params = {
            'db': 'pubmed',
            'id': ','.join(pmid_list),
            'retmode': 'xml'
        }
        
        url = self._build_url('efetch.fcgi', params)
        response = self._make_request(url)
        
        if not response:
            return []
        
        return self._parse_pubmed_xml(response)
    
    def _parse_pubmed_xml(self, xml_data: str) -> List[Dict]:
        """Parse PubMed XML response."""
        results = []
        
        try:
            root = ET.fromstring(xml_data)
            
            for article in root.findall('.//PubmedArticle'):
                record = {}
                
                # PMID
                pmid = article.find('.//PMID')
                if pmid is not None:
                    record['pmid'] = pmid.text
                
                # Title
                title = article.find('.//ArticleTitle')
                if title is not None:
                    record['title'] = title.text
                
                # Authors
                authors = []
                for author in article.findall('.//Author'):
                    last = author.find('.//LastName')
                    first = author.find('.//ForeName')
                    if last is not None:
                        name = last.text
                        if first is not None:
                            name = f"{first.text} {name}"
                        authors.append(name)
                record['authors'] = authors[:3]  # First 3 authors
                
                # Journal and date
                journal = article.find('.//Journal/Title')
                if journal is not None:
                    record['journal'] = journal.text
                
                pub_date = article.find('.//PubDate/Year')
                if pub_date is not None:
                    record['year'] = pub_date.text
                
                # Abstract
                abstract = article.find('.//AbstractText')
                if abstract is not None:
                    record['abstract'] = abstract.text
                
                results.append(record)
        
        except ET.ParseError as e:
            logger.error(f"XML Parse Error: {e}")
        
        return results
    
    def get_sra_from_pubmed(self, pmid: str) -> List[str]:
        """Get linked SRA accessions from a PubMed ID."""
        params = {
            'dbfrom': 'pubmed',
            'db': 'sra',
            'id': pmid,
            'retmode': 'json'
        }
        
        url = self._build_url('elink.fcgi', params)
        response = self._make_request(url)
        
        if not response:
            return []
        
        try:
            data = json.loads(response)
            linksets = data.get('linksets', [])
            
            sra_ids = []
            for linkset in linksets:
                for linksetdb in linkset.get('linksetdbs', []):
                    if linksetdb.get('dbto') == 'sra':
                        sra_ids.extend(linksetdb.get('links', []))
            
            return sra_ids
        
        except json.JSONDecodeError:
            return []
    
    def search_bioproject(self, query: str, retmax: int = 50) -> List[Dict]:
        """Search BioProject database for metagenomics projects."""
        params = {
            'db': 'bioproject',
            'term': query,
            'retmax': str(retmax),
            'retmode': 'json'
        }
        
        url = self._build_url('esearch.fcgi', params)
        logger.info(f"[BioProject Search] Query: {query}")
        
        response = self._make_request(url)
        if not response:
            return []
        
        try:
            data = json.loads(response)
            id_list = data.get('esearchresult', {}).get('idlist', [])
            count = data.get('esearchresult', {}).get('count', '0')
            logger.info(f"[BioProject] Found {count} projects, retrieving {len(id_list)}")
            
            if not id_list:
                return []
            
            return self._fetch_bioproject_details(id_list)
        
        except json.JSONDecodeError:
            logger.error("Error parsing BioProject results")
            return []
    
    def _fetch_bioproject_details(self, id_list: List[str]) -> List[Dict]:
        """Fetch BioProject details."""
        params = {
            'db': 'bioproject',
            'id': ','.join(id_list),
            'retmode': 'xml'
        }
        
        url = self._build_url('efetch.fcgi', params)
        response = self._make_request(url)
        
        if not response:
            return []
        
        return self._parse_bioproject_xml(response)
    
    def _parse_bioproject_xml(self, xml_data: str) -> List[Dict]:
        """Parse BioProject XML response."""
        results = []
        
        try:
            root = ET.fromstring(xml_data)
            
            for package in root.findall('.//Package'):
                record = {}
                
                # Project accession
                project = package.find('.//Project')
                if project is not None:
                    accn = project.find('.//ProjectID/ArchiveID')
                    if accn is not None:
                        record['accession'] = accn.get('accession', 'N/A')
                
                # Title and description
                descr = package.find('.//ProjectDescr')
                if descr is not None:
                    title = descr.find('.//Title')
                    if title is not None:
                        record['title'] = title.text
                    
                    desc = descr.find('.//Description')
                    if desc is not None:
                        record['description'] = desc.text
                
                # Project type
                proj_type = package.find('.//ProjectType')
                if proj_type is not None:
                    submission = proj_type.find('.//ProjectTypeSubmission')
                    if submission is not None:
                        record['project_type'] = submission.get('submission_type', 'N/A')
                
                # Organism
                organism = package.find('.//Organism')
                if organism is not None:
                    org_name = organism.find('.//OrganismName')
                    if org_name is not None:
                        record['organism'] = org_name.text
                
                results.append(record)
        
        except ET.ParseError as e:
            logger.error(f"XML Parse Error: {e}")
        
        return results
    
    def get_sra_from_bioproject(self, bioproject_acc: str) -> List[str]:
        """Get SRA run accessions associated with a BioProject."""
        # First, search SRA with the BioProject accession
        params = {
            'db': 'sra',
            'term': f'{bioproject_acc}[BioProject]',
            'retmax': '500',
            'retmode': 'json'
        }
        
        url = self._build_url('esearch.fcgi', params)
        response = self._make_request(url)
        
        if not response:
            return []
        
        try:
            data = json.loads(response)
            return data.get('esearchresult', {}).get('idlist', [])
        except json.JSONDecodeError:
            return []
    
    def build_sra_search_query(self, 
                              environment: Optional[str] = None,
                              pathogens: Optional[List[str]] = None,
                              host: Optional[str] = None,
                              has_short_reads: bool = True,
                              has_long_reads: bool = True) -> str:
        """Build SRA search query based on parameters."""
        query_parts = []
        
        # Base metagenomics filter
        query_parts.append("metagenome[All Fields]")
        
        # Environment filter
        if environment:
            query_parts.append(f'("{environment}"[Source] OR "{environment}"[All Fields])')
        
        # For metagenomes, pathogens are better searched in description/title
        if pathogens:
            pathogen_queries = []
            for p in pathogens:
                pathogen_queries.append(f'"{p}"[All Fields]')
            query_parts.append(f"({' OR '.join(pathogen_queries)})")
        
        # Host filter
        if host:
            query_parts.append(f'"{host}"[Organism]')
        
        # Platform filters
        platform_parts = []
        if has_short_reads:
            platform_parts.append('("ILLUMINA"[Platform] OR "BGISEQ"[Platform])')
        if has_long_reads:
            platform_parts.append('("OXFORD_NANOPORE"[Platform] OR "PACBIO_SMRT"[Platform])')
        
        # If both are requested, the intention is usually "either one" for standard search,
        # but if we want "hybrid" it's trickier in a single SRA query.
        # The current logic joins them with OR, which means "has short OR has long".
        # This matches the behavior of "I want samples that have sequencing data".

        # However, if 'hybrid_only' is used (which sets both has_short=True and has_long=True),
        # the user might expect "AND" logic if that were possible for a single run,
        # but for SRA, a "Run" is usually one platform.
        # So "hybrid" usually means "Same BioSample has runs from both platforms".
        #
        # For now, we will stick to the existing OR logic which retrieves candidates,
        # and we can't easily enforce "AND" at the search query level for *Runs*.
        #
        # BUT, wait! The previous implementation used OR.
        # If the user wants HYBRID, they probably want to ensure we search for BOTH types.
        #
        # Let's keep the OR logic here. If we wanted to enforce hybrid *strictly*,
        # we would need to post-process the results to ensure a BioSample has both.
        # Given the scope, enabling the flag to ensure both platforms are included in the OR query
        # is the first step.

        if platform_parts:
            query_parts.append(f"({' OR '.join(platform_parts)})")
        
        return " AND ".join(query_parts)
    
    def search_sra(self, query: str, retmax: int = 100, retstart: int = 0) -> Tuple[List[str], int]:
        """Search SRA database and return list of IDs and total count."""
        params = {
            'db': 'sra',
            'term': query,
            'retmax': str(retmax),
            'retstart': str(retstart),
            'retmode': 'json'
        }
        
        url = self._build_url('esearch.fcgi', params)
        if retstart == 0:
            logger.info(f"[SRA Search] Query: {query}")
        else:
            logger.info(f"[SRA Search] Fetching batch starting at {retstart}...")
        
        response = self._make_request(url)
        if not response:
            return [], 0
        
        try:
            data = json.loads(response)
            id_list = data.get('esearchresult', {}).get('idlist', [])
            count = int(data.get('esearchresult', {}).get('count', '0'))
            if retstart == 0:
                logger.info(f"[SRA] Found {count} total results, retrieving first batch")
            return id_list, count
        except json.JSONDecodeError:
            logger.error("Error parsing SRA search results")
            return [], 0
    
    def fetch_sra_details(self, id_list: List[str]) -> List[Dict]:
        """Fetch detailed information for given SRA IDs."""
        if not id_list:
            return []
        
        params = {
            'db': 'sra',
            'id': ','.join(id_list),
            'retmode': 'xml'
        }
        
        url = self._build_url('efetch.fcgi', params)
        response = self._make_request(url)
        
        if not response:
            return []
        
        return self._parse_sra_xml(response)
    
    def _parse_sra_xml(self, xml_data: str) -> List[Dict]:
        """Parse SRA XML response and extract relevant information."""
        results = []
        
        try:
            root = ET.fromstring(xml_data)
            
            for exp_pkg in root.findall('.//EXPERIMENT_PACKAGE'):
                record = {}
                
                # Accessions
                exp = exp_pkg.find('.//EXPERIMENT')
                if exp is not None:
                    record['experiment_accession'] = exp.get('accession', 'N/A')
                
                run = exp_pkg.find('.//RUN')
                if run is not None:
                    record['run_accession'] = run.get('accession', 'N/A')
                
                # Sample info
                sample = exp_pkg.find('.//SAMPLE')
                if sample is not None:
                    record['sample_accession'] = sample.get('accession', 'N/A')
                    sample_name = sample.find('.//SCIENTIFIC_NAME')
                    if sample_name is not None:
                        record['organism'] = sample_name.text
                
                # Study/BioProject
                study = exp_pkg.find('.//STUDY')
                if study is not None:
                    record['study_accession'] = study.get('accession', 'N/A')
                    ext_id = study.find('.//EXTERNAL_ID[@namespace="BioProject"]')
                    if ext_id is not None:
                        record['bioproject'] = ext_id.text
                
                # Platform
                platform = exp_pkg.find('.//PLATFORM')
                if platform is not None:
                    for child in platform:
                        record['platform'] = child.tag
                        instr = child.find('.//INSTRUMENT_MODEL')
                        if instr is not None:
                            record['instrument'] = instr.text
                
                # Library strategy
                lib_desc = exp_pkg.find('.//LIBRARY_DESCRIPTOR')
                if lib_desc is not None:
                    strategy = lib_desc.find('.//LIBRARY_STRATEGY')
                    if strategy is not None:
                        record['library_strategy'] = strategy.text
                
                # Sample attributes
                attributes = {}
                for attr in exp_pkg.findall('.//SAMPLE_ATTRIBUTE'):
                    tag = attr.find('.//TAG')
                    value = attr.find('.//VALUE')
                    if tag is not None and value is not None:
                        attributes[tag.text] = value.text
                
                record['sample_attributes'] = attributes
                
                if record:
                    results.append(record)
        
        except ET.ParseError as e:
            logger.error(f"XML Parse Error: {e}")
        
        return results
    
    def get_run_platforms_for_sample(self, sample_acc: str) -> List[str]:
        """Get all sequencing platforms used for a specific BioSample."""
        # Search SRA for runs associated with this sample
        params = {
            'db': 'sra',
            'term': f'{sample_acc}[Sample]',
            'retmax': '100',
            'retmode': 'json'
        }

        url = self._build_url('esearch.fcgi', params)
        response = self._make_request(url)

        if not response:
            return []

        try:
            data = json.loads(response)
            id_list = data.get('esearchresult', {}).get('idlist', [])
            if not id_list:
                return []

            # Fetch summary to get platforms (lighter than full XML)
            params = {
                'db': 'sra',
                'id': ','.join(id_list),
                'retmode': 'json',
                'version': '2.0'
            }
            url = self._build_url('esummary.fcgi', params)
            response = self._make_request(url)

            if not response:
                return []

            data = json.loads(response)
            result = data.get('result', {})
            platforms = set()

            for uid in id_list:
                if uid in result:
                    # Parse expxml or run info if available, but summary structure varies.
                    # In ESUMMARY 2.0 for SRA, platform is often in 'platform' field
                    item = result[uid]
                    if 'platform' in item and item['platform']:
                        platforms.add(item['platform'].upper())

                    # Sometimes it is in expxml
                    if 'expxml' in item:
                        # Simple regex to find platform - check for different capitalization
                        # <Platform instrument_model="...">ILLUMINA</Platform>
                        match = re.search(r'<Platform\s+instrument_model="[^"]+">([^<]+)</Platform>', item['expxml'], re.IGNORECASE)
                        if match:
                            platforms.add(match.group(1).upper())

                        # Fallback: sometimes it's just <Platform>NAME</Platform>
                        match_simple = re.search(r'<Platform>([^<]+)</Platform>', item['expxml'], re.IGNORECASE)
                        if match_simple:
                            platforms.add(match_simple.group(1).upper())

                    # Fallback to 'statistics' if present (sometimes has platform info?) - unlikely for summary

            return list(platforms)

        except json.JSONDecodeError:
            return []

    def validate_accession(self, accession: str) -> Tuple[bool, str]:
        """Validate if an accession exists and return its type."""
        # Determine database from accession prefix
        db_map = {
            'SRR': 'sra', 'ERR': 'sra', 'DRR': 'sra',
            'SRX': 'sra', 'ERX': 'sra', 'DRX': 'sra',
            'SAMN': 'biosample', 'SAME': 'biosample', 'SAMD': 'biosample',
            'PRJNA': 'bioproject', 'PRJEB': 'bioproject', 'PRJDB': 'bioproject'
        }
        
        prefix = accession[:3] if len(accession) >= 3 else None
        if accession.startswith('SAMN') or accession.startswith('SAME') or accession.startswith('SAMD'):
            prefix = accession[:4]
        elif accession.startswith('PRJNA') or accession.startswith('PRJEB') or accession.startswith('PRJDB'):
            prefix = accession[:5]
        
        db = db_map.get(prefix)
        if not db:
            return False, "Unknown accession format"
        
        params = {
            'db': db,
            'term': f'{accession}[Accession]',
            'retmode': 'json'
        }
        
        url = self._build_url('esearch.fcgi', params)
        response = self._make_request(url)
        
        if not response:
            return False, "API request failed"
        
        try:
            data = json.loads(response)
            count = int(data.get('esearchresult', {}).get('count', '0'))
            if count > 0:
                return True, f"Valid {db.upper()} accession"
            else:
                return False, f"Accession not found in {db.upper()}"
        except json.JSONDecodeError:
            return False, "Error parsing response"


def load_config(config_path: str) -> Dict:
    """Load configuration from YAML file."""
    try:
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        logger.error(f"Config file not found: {config_path}")
        sys.exit(1)
    except yaml.YAMLError as e:
        logger.error(f"Error parsing YAML: {e}")
        sys.exit(1)


def print_sra_results(results: List[Dict]):
    """Pretty print SRA search results."""
    if not results:
        print("No results found.")
        return
    
    print(f"\n{'='*80}")
    print(f"Found {len(results)} SRA samples:")
    print(f"{'='*80}\n")
    
    for i, record in enumerate(results, 1):
        print(f"--- Sample {i} ---")
        print(f"Run Accession:        {record.get('run_accession', 'N/A')}")
        print(f"Experiment Accession: {record.get('experiment_accession', 'N/A')}")
        print(f"Sample Accession:     {record.get('sample_accession', 'N/A')}")
        print(f"Study Accession:      {record.get('study_accession', 'N/A')}")
        if 'bioproject' in record:
            print(f"BioProject:           {record['bioproject']}")
        print(f"Organism:             {record.get('organism', 'N/A')}")
        print(f"Platform:             {record.get('platform', 'N/A')}")
        print(f"Instrument:           {record.get('instrument', 'N/A')}")
        print(f"Library Strategy:     {record.get('library_strategy', 'N/A')}")
        
        if record.get('sample_attributes'):
            print("\nSample Attributes:")
            for key, value in record['sample_attributes'].items():
                print(f"  {key}: {value}")
        
        print()


def print_bioproject_results(results: List[Dict]):
    """Pretty print BioProject results."""
    if not results:
        print("No BioProject results found.")
        return
    
    print(f"\n{'='*80}")
    print(f"Found {len(results)} BioProjects:")
    print(f"{'='*80}\n")
    
    for i, record in enumerate(results, 1):
        print(f"--- Project {i} ---")
        print(f"Accession:     {record.get('accession', 'N/A')}")
        print(f"Title:         {record.get('title', 'N/A')}")
        print(f"Type:          {record.get('project_type', 'N/A')}")
        print(f"Organism:      {record.get('organism', 'N/A')}")
        
        if 'description' in record and record['description']:
            desc = record['description'][:200] + "..." if len(record['description']) > 200 else record['description']
            print(f"Description:   {desc}")
        
        print()


def print_pubmed_results(results: List[Dict]):
    """Pretty print PubMed results."""
    if not results:
        print("No PubMed results found.")
        return
    
    print(f"\n{'='*80}")
    print(f"Found {len(results)} Publications:")
    print(f"{'='*80}\n")
    
    for i, record in enumerate(results, 1):
        print(f"--- Publication {i} ---")
        print(f"PMID:    {record.get('pmid', 'N/A')}")
        print(f"Title:   {record.get('title', 'N/A')}")
        
        if record.get('authors'):
            authors_str = ", ".join(record['authors'])
            if len(record['authors']) == 3:
                authors_str += ", et al."
            print(f"Authors: {authors_str}")
        
        if 'journal' in record:
            journal_str = record['journal']
            if 'year' in record:
                journal_str += f" ({record['year']})"
            print(f"Journal: {journal_str}")
        
        print(f"Link:    https://pubmed.ncbi.nlm.nih.gov/{record.get('pmid', '')}/")
        print()


def main():
    parser = argparse.ArgumentParser(
        description='Query NCBI for clinical metagenomics samples with short/long reads',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Search SRA directly
  %(prog)s --sra --environment blood --host "Homo sapiens" --max-results 10
  
  # Search BioProject for studies
  %(prog)s --bioproject --keywords "Klebsiella" "blood" "metagenome" --max-results 20
  
  # Search PubMed and find linked SRA data
  %(prog)s --pubmed --keywords "Klebsiella pneumoniae" "metagenomics" "bloodstream"
  
  # Get SRA runs from a specific BioProject
  %(prog)s --from-bioproject PRJNA12345
  
  # Validate accessions
  %(prog)s --validate SRR12345678 SAMN12345678 PRJNA12345
  
  # Use config file
  %(prog)s --config config.yaml --sra
        """
    )
    
    # Mode selection
    mode_group = parser.add_mutually_exclusive_group(required=False)
    mode_group.add_argument('--sra', action='store_true', 
                           help='Search SRA database (default)')
    mode_group.add_argument('--bioproject', action='store_true',
                           help='Search BioProject database')
    mode_group.add_argument('--pubmed', action='store_true',
                           help='Search PubMed for publications')
    mode_group.add_argument('--from-bioproject', metavar='PRJNA',
                           help='Get SRA runs from a BioProject accession')
    mode_group.add_argument('--from-pubmed', metavar='PMID',
                           help='Get SRA data linked to a PubMed ID')
    mode_group.add_argument('--validate', nargs='+', metavar='ACC',
                           help='Validate one or more accessions')
    
    # Search parameters
    parser.add_argument('--config', '-c', help='Path to YAML config file')
    parser.add_argument('--environment', '-e', 
                       help='Sample environment (e.g., blood, respiratory)')
    parser.add_argument('--pathogens', '-p', nargs='+', 
                       help='List of pathogen names to search for')
    parser.add_argument('--host', '-H', 
                       help='Host organism (e.g., "Homo sapiens")')
    parser.add_argument('--keywords', '-k', nargs='+',
                       help='Keywords for BioProject/PubMed search')
    
    # API parameters
    parser.add_argument('--email', default='user@example.com', 
                       help='Email for NCBI (required by API)')
    parser.add_argument('--api-key', help='NCBI API key for higher rate limits')
    
    # Output parameters
    parser.add_argument('--max-results', '-m', type=int, default=20, 
                       help='Maximum results to retrieve')
    parser.add_argument('--no-short-reads', action='store_true', 
                       help='Exclude short read requirement')
    parser.add_argument('--no-long-reads', action='store_true', 
                       help='Exclude long read requirement')
    parser.add_argument('--hybrid-only', action='store_true',
                       help='Require both short and long reads for the same experiment (simulated by requiring both)')
    parser.add_argument('--output', '-o', help='Output file (JSON format)')
    parser.add_argument('--log', help='Log file path')
    parser.add_argument('--get-sra', action='store_true',
                       help='For BioProject/PubMed: also fetch linked SRA data')
    
    args = parser.parse_args()

    # Setup logging
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    if args.log:
        logging.basicConfig(filename=args.log, level=logging.INFO, format=log_format)
        # Also log to stderr
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        formatter = logging.Formatter(log_format)
        console.setFormatter(formatter)
        logger.addHandler(console)
    else:
        logging.basicConfig(level=logging.INFO, format='%(message)s')
    
    # Default to SRA search if no mode specified
    if not any([args.sra, args.bioproject, args.pubmed, args.from_bioproject, 
                args.from_pubmed, args.validate]):
        args.sra = True
    
    # Load config if provided
    config = {}
    if args.config:
        config = load_config(args.config)
    
    # Get parameters (CLI overrides config)
    environment = args.environment or config.get('environment')
    pathogens = args.pathogens or config.get('pathogens')
    host = args.host or config.get('host')
    keywords = args.keywords or config.get('keywords', [])
    email = args.email or config.get('email', 'user@example.com')
    api_key = args.api_key or config.get('api_key')
    
    has_short = not args.no_short_reads
    has_long = not args.no_long_reads
    
    if args.hybrid_only:
        has_short = True
        has_long = True

    # Initialize tool
    tool = EntrezQueryTool(email=email, api_key=api_key)
    
    results = []
    
    # Validation mode
    if args.validate:
        print(f"\n{'='*80}")
        print("Validating Accessions:")
        print(f"{'='*80}\n")
        
        for acc in args.validate:
            is_valid, message = tool.validate_accession(acc)
            status = "✓ VALID" if is_valid else "✗ INVALID"
            print(f"{acc}: {status} - {message}")
        return
    
    # From BioProject mode
    if args.from_bioproject:
        logger.info(f"Fetching SRA runs from BioProject: {args.from_bioproject}")
        sra_ids = tool.get_sra_from_bioproject(args.from_bioproject)
        
        if sra_ids:
            logger.info(f"Found {len(sra_ids)} associated SRA runs")
            results = tool.fetch_sra_details(sra_ids[:args.max_results])
            print_sra_results(results)
        else:
            logger.info("No SRA data found for this BioProject")
        
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Results saved to {args.output}")
        return
    
    # From PubMed mode
    if args.from_pubmed:
        logger.info(f"Fetching SRA data linked to PMID: {args.from_pubmed}")
        sra_ids = tool.get_sra_from_pubmed(args.from_pubmed)
        
        if sra_ids:
            logger.info(f"Found {len(sra_ids)} linked SRA runs")
            results = tool.fetch_sra_details(sra_ids[:args.max_results])
            print_sra_results(results)
        else:
            logger.info("No SRA data linked to this publication")
        
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Results saved to {args.output}")
        return
    
    # PubMed search mode
    if args.pubmed:
        if not keywords:
            logger.error("Error: --keywords required for PubMed search")
            sys.exit(1)
        
        query = " AND ".join([f'"{k}"' for k in keywords])
        results = tool.search_pubmed(query, retmax=args.max_results)
        print_pubmed_results(results)
        
        if args.get_sra and results:
            logger.info("Fetching linked SRA data for publications...")
            all_sra = []
            for pub in results[:5]:  # Limit to first 5 to avoid too many requests
                pmid = pub.get('pmid')
                if pmid:
                    sra_ids = tool.get_sra_from_pubmed(pmid)
                    if sra_ids:
                        all_sra.extend(sra_ids)

            if all_sra:
                # Remove duplicates
                all_sra = list(set(all_sra))
                logger.info(f"Found {len(all_sra)} unique SRA runs linked to publications")
                sra_results = tool.fetch_sra_details(all_sra[:args.max_results])
                print_sra_results(sra_results)
            else:
                logger.info("No linked SRA data found")

        if args.output:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Results saved to {args.output}")
        return

    # SRA search mode (default)
    if args.sra:
        query = tool.build_sra_search_query(
            environment=environment,
            pathogens=pathogens,
            host=host,
            has_short_reads=has_short,
            has_long_reads=has_long
        )

        final_details = []

        if args.hybrid_only:
            logger.info("Filtering for hybrid samples (both Short and Long reads)...")
            processed_samples = set()
            valid_samples = set()

            batch_size = 50
            start = 0
            max_search_limit = 1000  # Safety limit to prevent infinite loops

            while len(valid_samples) < args.max_results and start < max_search_limit:
                # Fetch batch of IDs
                results, total_count = tool.search_sra(query, retmax=batch_size, retstart=start)

                if not results:
                    break

                # Fetch details for this batch
                batch_details = tool.fetch_sra_details(results)

                for record in batch_details:
                    # If we have enough, stop processing
                    if len(valid_samples) >= args.max_results:
                        break

                    sample_acc = record.get('sample_accession')
                    if not sample_acc or sample_acc == 'N/A':
                        continue

                    # Check each sample only once
                    if sample_acc in processed_samples:
                        if sample_acc in valid_samples:
                            final_details.append(record)
                        continue

                    processed_samples.add(sample_acc)
                    logger.info(f"Checking sample {sample_acc}...")

                    platforms = tool.get_run_platforms_for_sample(sample_acc)

                    has_illumina_bgi = any(p in ['ILLUMINA', 'BGISEQ'] for p in platforms)
                    has_long = any(p in ['OXFORD_NANOPORE', 'PACBIO_SMRT'] for p in platforms)

                    if has_illumina_bgi and has_long:
                        logger.info(f"Sample {sample_acc}: HYBRID FOUND!")
                        valid_samples.add(sample_acc)
                        final_details.append(record)
                    else:
                        logger.info(f"Sample {sample_acc}: Platforms: {', '.join(platforms)}")

                start += batch_size
                if start >= total_count:
                    break

            logger.info(f"Found {len(valid_samples)} hybrid samples after checking {len(processed_samples)} candidates.")

        else:
            # Standard search
            results, _ = tool.search_sra(query, retmax=args.max_results)
            final_details = tool.fetch_sra_details(results)

        print_sra_results(final_details)

        if args.output:
            with open(args.output, 'w') as f:
                json.dump(final_details, f, indent=2)
            logger.info(f"Results saved to {args.output}")


if __name__ == "__main__":
    main()
