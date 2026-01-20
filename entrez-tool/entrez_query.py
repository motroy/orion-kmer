#!/usr/bin/env python3
"""
NCBI Entrez CLI tool for querying samples with both short and long read data available.

Enhanced with BioProject search, PubMed links, and validation mode.
Refactored to use pysradb for SRA metadata and ENA links.
Refactored to use metapub for PubMed queries.
"""

import argparse
import sys
import yaml
import logging
from typing import List, Optional, Dict, Tuple, Set
from urllib.parse import quote
from urllib.request import urlopen
from urllib.error import HTTPError, URLError
import xml.etree.ElementTree as ET
import time
import json
import re
import pandas as pd
from pysradb.sraweb import SRAweb
from metapub import PubMedFetcher

# Configure logging
logger = logging.getLogger('entrez_tool')
logger.setLevel(logging.INFO)

class EntrezQueryTool:
    """Tool for querying NCBI Entrez databases for samples."""
    
    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    def __init__(self, email: str = "user@example.com", api_key: Optional[str] = None):
        self.email = email
        self.api_key = api_key
        self.delay = 0.34 if not api_key else 0.1
        self.sra_web = SRAweb()
        self.pubmed_fetcher = PubMedFetcher(apikey=api_key) if api_key else PubMedFetcher()
    
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
        """Search PubMed for publications using metapub."""
        logger.info(f"[PubMed Search] Query: {query}")
        try:
            pmids = self.pubmed_fetcher.pmids_for_query(query, retmax=retmax)
            logger.info(f"[PubMed] Found {len(pmids)} publications")
            
            if not pmids:
                return []
            
            results = []
            for pmid in pmids:
                try:
                    article = self.pubmed_fetcher.article_by_pmid(pmid)
                    if article:
                        results.append(article)
                except Exception as e:
                    logger.warning(f"Error fetching article for PMID {pmid}: {e}")
            return results
        except Exception as e:
            logger.error(f"PubMed search failed: {e}")
            return []
    
    def get_sra_from_pubmed(self, pmid: str) -> List[str]:
        """Get linked SRA accessions from a PubMed ID."""
        # metapub doesn't seem to have direct elink support for SRA, so we keep using eutils directly.
        params = {
            'dbfrom': 'pubmed',
            'db': 'sra',
            'id': pmid,
            'retmode': 'json'
        }
        url = self._build_url('elink.fcgi', params)
        response = self._make_request(url)
        if not response: return []
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
        if not response: return []
        try:
            data = json.loads(response)
            id_list = data.get('esearchresult', {}).get('idlist', [])
            count = data.get('esearchresult', {}).get('count', '0')
            logger.info(f"[BioProject] Found {count} projects, retrieving {len(id_list)}")
            if not id_list: return []
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
        if not response: return []
        return self._parse_bioproject_xml(response)
    
    def _parse_bioproject_xml(self, xml_data: str) -> List[Dict]:
        """Parse BioProject XML response."""
        results = []
        try:
            root = ET.fromstring(xml_data)
            for package in root.findall('.//Package'):
                record = {}
                project = package.find('.//Project')
                if project is not None:
                    accn = project.find('.//ProjectID/ArchiveID')
                    if accn is not None: record['accession'] = accn.get('accession', 'N/A')
                descr = package.find('.//ProjectDescr')
                if descr is not None:
                    title = descr.find('.//Title')
                    if title is not None: record['title'] = title.text
                    desc = descr.find('.//Description')
                    if desc is not None: record['description'] = desc.text
                proj_type = package.find('.//ProjectType')
                if proj_type is not None:
                    submission = proj_type.find('.//ProjectTypeSubmission')
                    if submission is not None: record['project_type'] = submission.get('submission_type', 'N/A')
                organism = package.find('.//Organism')
                if organism is not None:
                    org_name = organism.find('.//OrganismName')
                    if org_name is not None: record['organism'] = org_name.text
                results.append(record)
        except ET.ParseError as e:
            logger.error(f"XML Parse Error: {e}")
        return results
    
    def get_sra_from_bioproject(self, bioproject_acc: str) -> List[str]:
        """Get SRA run accessions associated with a BioProject."""
        params = {
            'db': 'sra',
            'term': f'{bioproject_acc}[BioProject]',
            'retmax': '500',
            'retmode': 'json'
        }
        url = self._build_url('esearch.fcgi', params)
        response = self._make_request(url)
        if not response: return []
        try:
            data = json.loads(response)
            return data.get('esearchresult', {}).get('idlist', [])
        except json.JSONDecodeError:
            return []
    
    # --- SRA Search Logic using pysradb ---

    def build_sra_search_query(self, 
                              environment: Optional[str] = None,
                              pathogens: Optional[List[str]] = None,
                              host: Optional[str] = None,
                              keywords: Optional[List[str]] = None,
                              has_short_reads: bool = True,
                              has_long_reads: bool = True) -> str:
        """Build SRA search query based on parameters."""
        query_parts = []
        
        # Keywords filter
        if keywords:
            keyword_queries = []
            for k in keywords:
                keyword_queries.append(f'"{k}"[All Fields]')
            query_parts.append(f"({' AND '.join(keyword_queries)})")
        
        # Environment filter
        if environment:
            query_parts.append(f'("{environment}"[Source] OR "{environment}"[All Fields])')
        
        # Pathogens
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
        
        if platform_parts:
            query_parts.append(f"({' OR '.join(platform_parts)})")
        
        return " AND ".join(query_parts)
    
    def search_sra(self, query: str, retmax: int = 100, retstart: int = 0) -> Tuple[List[str], int]:
        """Search SRA database and return list of UIDs and total count using Entrez API."""
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

    def get_accessions_from_uids(self, uids: List[str]) -> List[str]:
        """Convert SRA UIDs to Accessions (SRR/ERR/DRR) using eSummary."""
        if not uids:
            return []
        
        params = {
            'db': 'sra',
            'id': ','.join(uids),
            'retmode': 'json'
        }
        
        url = self._build_url('esummary.fcgi', params)
        response = self._make_request(url)
        if not response:
            return []
        
        accessions = []
        try:
            data = json.loads(response)
            result = data.get('result', {})
            for uid in uids:
                if uid in result:
                    item = result[uid]
                    # Extract Run Accession from 'runs' string if available
                    runs_str = item.get('runs', '')
                    matches = re.findall(r'acc="([SED]RR\d+)"', runs_str)
                    if matches:
                        accessions.extend(matches)
                    else:
                        expxml = item.get('expxml', '')
                        match_exp = re.search(r'Experiment\s+acc="([SED]RX\d+)"', expxml, re.IGNORECASE)
                        if match_exp:
                            accessions.append(match_exp.group(1))
        except Exception as e:
            logger.error(f"Error converting UIDs to Accessions: {e}")
        
        return list(set(accessions))

    def fetch_sra_details(self, id_list: List[str]) -> List[Dict]:
        """Fetch detailed information for given SRA UIDs using pysradb."""
        if not id_list:
            return []

        # 1. Convert UIDs to Accessions
        accessions = self.get_accessions_from_uids(id_list)
        if not accessions:
            return []

        logger.info(f"Fetching metadata for {len(accessions)} accessions using pysradb...")

        # 2. Fetch metadata using pysradb
        try:
            df = self.sra_web.sra_metadata(accessions, detailed=True)
        except Exception as e:
            logger.error(f"pysradb metadata fetch failed: {e}")
            return []

        # 3. Convert DataFrame to List[Dict]
        return self._dataframe_to_results(df)

    def _dataframe_to_results(self, df: pd.DataFrame) -> List[Dict]:
        """Convert pysradb DataFrame to list of dictionaries matching legacy output."""
        results = []
        if df is None or df.empty:
            return results

        for _, row in df.iterrows():
            record = {}
            # Map columns
            record['run_accession'] = row.get('run_accession', 'N/A')
            record['experiment_accession'] = row.get('experiment_accession', 'N/A')
            record['sample_accession'] = row.get('sample_accession', 'N/A')
            record['study_accession'] = row.get('study_accession', 'N/A')
            record['bioproject'] = row.get('bioproject', 'N/A')
            record['organism'] = row.get('organism_name', 'N/A')

            # Platform/Instrument
            record['instrument'] = row.get('instrument_model', 'N/A')
            record['platform'] = "N/A"
            if record['instrument']:
                instr_str = str(record['instrument'])
                if 'Illumina' in instr_str:
                    record['platform'] = 'ILLUMINA'
                elif 'MinION' in instr_str or 'PromethION' in instr_str or 'Nanopore' in instr_str:
                    record['platform'] = 'OXFORD_NANOPORE'
                elif 'PacBio' in instr_str or 'Sequel' in instr_str:
                    record['platform'] = 'PACBIO_SMRT'

            record['library_strategy'] = row.get('library_strategy', 'N/A')

            # ENA Links
            ena_links = []
            for col in ['ena_fastq_ftp', 'ena_fastq_ftp_1', 'ena_fastq_ftp_2']:
                if col in row and not pd.isna(row[col]):
                    ena_links.append(row[col])

            if ena_links:
                record['ena_links'] = ena_links

            results.append(record)
        return results

    def get_run_platforms_for_sample(self, sample_acc: str) -> List[str]:
        """Get all sequencing platforms used for a specific BioSample using pysradb."""
        try:
            df = self.sra_web.sra_metadata(sample_acc, detailed=False)
            if df is None or df.empty:
                return []

            platforms = set()
            if 'instrument_model' in df.columns:
                for instr in df['instrument_model']:
                    if pd.isna(instr): continue
                    instr_upper = str(instr).upper()
                    if 'ILLUMINA' in instr_upper: platforms.add('ILLUMINA')
                    elif 'BGI' in instr_upper: platforms.add('BGISEQ')
                    elif 'NANOPORE' in instr_upper or 'MINION' in instr_upper or 'PROMETHION' in instr_upper: platforms.add('OXFORD_NANOPORE')
                    elif 'PACBIO' in instr_upper or 'SEQUEL' in instr_upper: platforms.add('PACBIO_SMRT')
                    else: platforms.add(instr_upper) # Fallback
            return list(platforms)
        except Exception as e:
            logger.error(f"Error fetching platforms for sample {sample_acc}: {e}")
            return []

    def validate_accession(self, accession: str) -> Tuple[bool, str]:
        """Validate if an accession exists and return its type."""
        # Use existing logic as it's simple
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
        
        if record.get('ena_links'):
            print("ENA Links:")
            for link in record['ena_links']:
                print(f"  {link}")

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


def print_pubmed_results(results: List[any]):
    """Pretty print PubMed results (metapub objects)."""
    if not results:
        print("No PubMed results found.")
        return
    
    print(f"\n{'='*80}")
    print(f"Found {len(results)} Publications:")
    print(f"{'='*80}\n")
    
    for i, article in enumerate(results, 1):
        print(f"--- Publication {i} ---")
        print(f"PMID:    {article.pmid}")
        print(f"Title:   {article.title}")
        
        if article.authors:
            authors_str = ", ".join(article.authors)
            if len(article.authors) > 5:
                authors_str = ", ".join(article.authors[:5]) + ", et al."
            print(f"Authors: {authors_str}")
        
        journal_str = article.journal or ""
        if article.year:
            journal_str += f" ({article.year})"
        print(f"Journal: {journal_str}")
        
        print(f"Link:    https://pubmed.ncbi.nlm.nih.gov/{article.pmid}/")
        print()


def main():
    parser = argparse.ArgumentParser(
        description='Query NCBI for samples with short/long reads',
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
                       help='Keywords for search')
    
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
        # Optimization: When looking for hybrid, search for the rarer technology (Long Reads)
        # to filter candidates, then check for Short Reads.
        has_short = False
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
            for article in results[:5]:  # Limit to first 5
                pmid = article.pmid
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
            # Need to serialize metapub objects manually or Convert them
            serializable_results = []
            for article in results:
                serializable_results.append({
                    'pmid': article.pmid,
                    'title': article.title,
                    'authors': article.authors,
                    'journal': article.journal,
                    'year': article.year,
                    'abstract': article.abstract
                })
            with open(args.output, 'w') as f:
                json.dump(serializable_results, f, indent=2)
            logger.info(f"Results saved to {args.output}")
        return

    # SRA search mode (default)
    if args.sra:
        query = tool.build_sra_search_query(
            environment=environment,
            pathogens=pathogens,
            host=host,
            keywords=keywords,
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
            max_search_limit = 1000  # Safety limit

            while len(valid_samples) < args.max_results and start < max_search_limit:
                # Fetch batch of UIDs
                results_uids, total_count = tool.search_sra(query, retmax=batch_size, retstart=start)

                if not results_uids:
                    break

                # Convert to accessions (needed for pysradb)
                # Note: this conversion might return SRR or SRX.
                # get_run_platforms_for_sample expects Sample Accession usually, but can work with others?
                # Actually, we need to know the Sample Accession for each Run to check for hybrid.

                # Let's get metadata for these runs using pysradb first
                batch_details = tool.fetch_sra_details(results_uids)

                for record in batch_details:
                    if len(valid_samples) >= args.max_results:
                        break

                    sample_acc = record.get('sample_accession')
                    if not sample_acc or sample_acc == 'N/A':
                        continue

                    # Check each sample only once
                    if sample_acc in processed_samples:
                        # If already valid, we might want to add this run too?
                        # The original logic just added the record if valid.
                        if sample_acc in valid_samples:
                            final_details.append(record)
                        continue

                    processed_samples.add(sample_acc)
                    logger.info(f"Checking sample {sample_acc}...")

                    # Use pysradb to check platforms for this sample
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
            results_uids, _ = tool.search_sra(query, retmax=args.max_results)
            final_details = tool.fetch_sra_details(results_uids)

        print_sra_results(final_details)

        if args.output:
            with open(args.output, 'w') as f:
                json.dump(final_details, f, indent=2)
            logger.info(f"Results saved to {args.output}")


if __name__ == "__main__":
    main()
