import datetime
import os
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Set, Tuple

import requests
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from bs4 import BeautifulSoup
from dotenv import load_dotenv

from bio_files_processor import FastaRecord


load_dotenv("tg_api_token.env")
TELEGRAM_API_URL = f"https://api.telegram.org/bot{os.getenv('TG_API_TOKEN')}/sendMessage"


def filter_fastq(input_path: str, 
                 output_filename: str = None, 
                 gc_bounds: Tuple[int, int] = (0, 100), 
                 length_bounds: Tuple[int, int] = (0, 2**32), 
                 quality_threshold: int = 0) -> None:
    '''
    Filters sequences in a FASTQ file based on specified criteria and saves the filtered sequences to a new FASTQ file.

    Parameters:
        input_path (str): The path to the input FASTQ file.
        output_filename (str, optional): The name of the output FASTQ file. If not provided, 
                                         a default name will be generated based on the input file name.
        gc_bounds (Tuple[int, int], optional): A tuple specifying the lower and upper bounds for GC content percentage.
                                        Default is (0, 100).
        length_bounds (Tuple[int, int], optional): A tuple specifying the lower and upper bounds for sequence length.
                                        Default is (0, 2**32).
        quality_threshold (int, optional): The minimum average quality score threshold for a sequence to be included.
                                        Default is 0.
    Returns:
        None: The filtered sequences are saved to a new FASTQ file.
    '''
    with open(input_path, "r") as handle:
        records = list(SeqIO.parse(handle, "fastq"))

    filtered_records = []
    for record in records:
        gc_content = gc_fraction(record.seq) * 100
        seq_length = len(record.seq)
        avg_quality = sum(record.letter_annotations["phred_quality"]) / seq_length

        if (gc_bounds[0] <= gc_content <= gc_bounds[1] and
            length_bounds[0] <= seq_length <= length_bounds[1] and
            avg_quality >= quality_threshold):
            filtered_records.append(record)

    if not os.path.exists("fastq_filtrator_results"):
        os.makedirs("fastq_filtrator_results")

    if output_filename is None:
        output_filename = input_path.split("/")[-1].split(".")[0] + "_filtered.fastq"

    output_path = "fastq_filtrator_results/" + output_filename
    SeqIO.write(filtered_records, output_path, "fastq")

    print(f"Filtered sequences saved to {output_path}")



class BiologicalSequence(ABC):
    """
    Abstract base class representing a biological sequence.
    """
    ALPHABET: Set[str] = None

    def __init__(self, seq: str = None) -> None:
        self.seq = seq

    @abstractmethod
    def __len__(self) -> int:
        pass

    @abstractmethod
    def __getitem__(self, index) -> str:
        pass

    @abstractmethod
    def __repr__(self) -> str:
        pass

    @abstractmethod
    def check_alphabet(self) -> bool:
        pass


class InvalidAlphabetError(Exception):
    """Exception raised when the sequence contains characters outside the expected alphabet."""
    pass


class NucleicAcidSequence(BiologicalSequence):
    """
    Represents a nucleic acid sequence.
    """
    def __len__(self) -> int:
        return len(self.seq)

    def __getitem__(self, index) -> str:
        return self.seq[index]

    def __repr__(self) -> str:
        return self.seq

    def check_alphabet(self) -> bool:
        """
        Check if the sequence contains valid nucleotide characters.
        """
        if not set(self.seq.upper()).issubset(self.ALPHABET):
            raise InvalidAlphabetError("Invalid characters found in the sequence.")
        return True
    
    def complement(self) -> 'NucleicAcidSequence':
        """
        Get the complement of the sequence.
        """
        comp_seq = self.seq.translate(str.maketrans(self.COMPLEMENT_MAP))
        return type(self)(comp_seq)
    
    def reverse_complement(self) -> 'NucleicAcidSequence':
        """
        Get the reverse complement of the sequence.
        """
        reverse_seq = self.complement().seq[::-1]
        return type(self)(reverse_seq)
    
    def gc_content(self) -> float:
        """
        Calculate the GC content of the sequence.
        """
        gc_count = self.seq.count('G') + self.seq.count('C')
        return gc_count / len(self.seq) * 100
    

class DNASequence(NucleicAcidSequence):
    """
    Represents a DNA sequence.
    """
    ALPHABET: Set[str] = {'A', 'T', 'G', 'C'}
    COMPLEMENT_MAP: dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                            'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    
    def transcribe(self) -> 'RNASequence':
        """
        Transcribe the DNA sequence to RNA.
        """
        return RNASequence(self.seq.replace('T', 'U'))


class RNASequence(NucleicAcidSequence):
    """
    Represents an RNA sequence.
    """
    ALPHABET: Set[str] = {'A', 'G', 'C', 'U'}
    COMPLEMENT_MAP: dict = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
                            'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}


class AminoAcidSequence(BiologicalSequence):
    """
    Represents an amino acid sequence.
    """
    ALPHABET: Set[str] = {"A", "C", "D", "E", "F", "G", "H", "I","K", "L", 
                          "M", "N","P", "Q", "R", "S", "T", "V", "W", "Y"}
    HYDROPHOBIC_AMINOACIDS: Set[str] = {"A", "V", "L", "I", "M", "F", "Y", "W"}
    
    def __len__(self) -> int:
        return len(self.seq)

    def __getitem__(self, index) -> str:
        return self.seq[index]

    def __repr__(self) -> str:
        return self.seq

    def check_alphabet(self) -> bool:
        """
        Check if the sequence contains valid amino acid characters.
        """
        if not set(self.seq.upper()).issubset(self.ALPHABET):
            raise InvalidAlphabetError("Invalid characters found in the sequence.")
        return True
    
    def compute_hydrophobicity(self) -> float:
        """
        Compute the hydrophobicity of the sequence.
        """
        hydrophobic_count = sum(1 for aa in self.seq if aa.upper() in self.HYDROPHOBIC_AMINOACIDS)
        return (hydrophobic_count / len(self.seq)) * 100


def format_execution_time(start_time, end_time):
    """
    Format the execution time between two datetime objects into a human-readable string.

    Parameters:
        start_time (datetime.datetime): The start time of execution.
        end_time (datetime.datetime): The end time of execution.

    Returns:
        str: A formatted string representing the execution time.
    """
    execution_time = end_time - start_time
    if execution_time.days >= 1:
        days = execution_time.days
        hours, remainder = divmod(execution_time.seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        return f"{days} days, {hours:02}:{minutes:02}:{seconds:02}"
    else:
        hours, remainder = divmod(execution_time.seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        return f"{hours:02}:{minutes:02}:{seconds:02}"


def send_telegram_message(chat_id, message, url):
    """
    Send a message to a Telegram chat using the Telegram Bot API.

    Parameters:
        chat_id (str): The ID of the Telegram chat.
        message (str): The message to be sent.
        url (str): The URL of the Telegram Bot API.

    Returns:
        None
    """
    data = {
        "chat_id": chat_id,
        "text": message
    }
    response = requests.post(url, data=data)
    if response.status_code != 200:
        print("Failed to send Telegram message.")
        print(response.text)


def telegram_logger(chat_id):
    """
    Decorator function to log the execution of a function and send a Telegram message upon success or failure.

    Parameters:
        chat_id (str): The ID of the Telegram chat where the messages will be sent.

    Returns:
        function: A decorator function.
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            start_time = datetime.datetime.now()

            try:
                result = func(*args, **kwargs)
                end_time = datetime.datetime.now()
                execution_time = format_execution_time(start_time, end_time)
                message = f"💅Function {func.__name__} finished successfully.\n"
                message += f"⏳Execution time: {execution_time}\n"
                send_telegram_message(chat_id, message, TELEGRAM_API_URL)

                return result
            
            except Exception as e:
                end_time = datetime.datetime.now()
                execution_time = format_execution_time(start_time, end_time)
                message = f"😞Function {func.__name__} failed with an exception:\n"
                message += f"{type(e).__name__}: {str(e)}\n"
                send_telegram_message(chat_id, message, TELEGRAM_API_URL)
                raise

        return wrapper

    return decorator


class Exon:
    """Represents an exon in a nucleotide sequence.

    Attributes:
        number (str): The exon number.
        type (str): The type of exon:
            - Init = Initial exon (ATG to 5' splice site)
            - Intr = Internal exon (3' splice site to 5' splice site)
            - Term = Terminal exon (3' splice site to stop codon)
            - Sngl = Single-exon gene (ATG to stop)
            - Prom = Promoter (TATA box / initation site)
            - PlyA = poly-A signal (consensus: AATAAA)
        start (int): The starting position of the exon.
        end (int): The ending position of the exon.
    """
    def __init__(self, number, exon_type, start, end):
        self.number = number
        self.type = exon_type
        self.start = start
        self.end = end

    def __repr__(self) -> str:
        return f"Exon {self.number}: Type={self.type}, Start={self.start}, End={self.end}"
    

class Intron:
    """Represents an intron in a nucleotide sequence.

    Attributes:
        number (int): The intron number.
        start (int): The starting position of the intron.
        end (int): The ending position of the intron.
    """
    def __init__(self, number, start, end):
        self.number = number
        self.start = start
        self.end = end

    def __repr__(self) -> str:
        return f"Intron {self.number}: Start={self.start}, End={self.end}"


@dataclass
class GenscanOutput:
    """Represents the output of a GENSCAN prediction.

    Attributes:
        status (int): The status code of the prediction response.
        cds_list (List[ProteinSequence]): List of predicted protein sequences.
        intron_list (List[Intron]): List of predicted introns.
        exon_list (List[Exon]): List of predicted exons.
    """
    status: int
    cds_list: list
    intron_list: list
    exon_list: list


def extract_introns(exons, sequence_length):
    """Extracts introns from a list of exons and sequence length.

    Args:
        exons (List[Exon]): List of exon objects.
        sequence_length (int): Length of the nucleotide sequence.

    Returns:
        List[Intron]: List of intron objects.
    """
    introns = []
    for i in range(len(exons) - 1):
        intron_start = exons[i].end + 1
        intron_end = exons[i + 1].start - 1
        intron = Intron(i + 1, intron_start, intron_end)
        introns.append(intron)
    return introns


def run_genscan(sequence: str=None, 
                sequence_file_path: str=None, 
                organism: str="Vertebrate", 
                exon_cutoff: float=1.00, 
                sequence_name: str=""):
    """Runs GENSCAN prediction and extracts the output.

    Args:
        sequence (str): Nucleotide sequence.
        sequence_file_path (str): Path to the file containing the nucleotide sequence.
        organism (str): Organism type for prediction.
        exon_cutoff (float): Exon cutoff value.
        sequence_name (str): Name of the sequence.

    Returns:
        GenscanOutput: Object containing prediction output.
    """
    site_url = "http://hollywood.mit.edu/cgi-bin/genscanw_py.cgi"
    form_data = {
    "-o": organism,
    "-n": sequence_name,
    "-u": open(sequence_file_path, 'rb'),
    "-e": str(exon_cutoff).encode(),
    "-s": sequence,
    "-p": "Predicted peptides only"
    }
    response = requests.post(site_url, files=form_data)
    status = response.status_code
    soup = BeautifulSoup(response.content, 'html.parser')
    content = soup.find('pre').text.split('\n')  
    cds_list = []
    exons = []
    current_sequence = ""
    current_name = None

    for line in content:
        line = line.strip()
        if line.startswith('Sequence'):
            sequence_length = int(line.split(' : ')[1].replace(' bp',''))
        elif line.startswith('1'):  
            parts = line.split()
            number = parts[0].split('.')[1].replace('0','')
            exon_type = parts[1]
            start = int(parts[3])
            end = int(parts[4])
            exon = Exon(number, exon_type, start, end)
            exons.append(exon)
        elif line.startswith('>'):  
            if current_name and current_sequence:
                sequence = FastaRecord(id=current_name, seq=current_sequence, description='')
                cds_list.append(sequence)
            current_name = line[1:]
            current_sequence = ""
        else:
            current_sequence += line.strip()
    
    if current_name and current_sequence:
        sequence = FastaRecord(id=current_name, seq=current_sequence, description='')
        cds_list.append(sequence)

    introns = extract_introns(exons, sequence_length)

    return GenscanOutput(status, cds_list=cds_list, intron_list=introns, exon_list=exons)
