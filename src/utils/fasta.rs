use bio::io::fasta;
use std::collections::HashMap;
use std::path::Path;

pub(crate) fn read_reference<P: AsRef<Path> + std::fmt::Debug>(
    fasta_path: P,
) -> HashMap<String, Vec<u8>> {
    let reader = fasta::Reader::from_file(fasta_path).expect("Error reading FASTA file");
    let mut genome_map = HashMap::new();

    for result in reader.records() {
        let record = result.expect("Error reading record");
        let id = record.id().to_string();
        let seq = record.seq().to_vec();
        genome_map.insert(id, seq);
    }

    genome_map
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_temp_fasta_file(contents: &str) -> NamedTempFile {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "{}", contents).unwrap();
        temp_file
    }

    #[test]
    fn reads_single_record_correctly() {
        let fasta_contents = ">seq1\nACGT";
        let temp_file = create_temp_fasta_file(fasta_contents);
        let result = read_reference(temp_file.path());
        assert_eq!(result.len(), 1);
        assert_eq!(result["seq1"], b"ACGT".to_vec());
    }

    #[test]
    fn reads_multiple_records_correctly() {
        let fasta_contents = ">seq1\nACGT\n>seq2\nTGCA";
        let temp_file = create_temp_fasta_file(fasta_contents);
        let result = read_reference(temp_file.path());
        assert_eq!(result.len(), 2);
        assert_eq!(result["seq1"], b"ACGT".to_vec());
        assert_eq!(result["seq2"], b"TGCA".to_vec());
    }

    #[test]
    #[should_panic(expected = "Error reading FASTA file")]
    fn panics_on_nonexistent_file() {
        let fasta_path = "nonexistent.fasta";
        let _ = read_reference(fasta_path);
    }

    #[test]
    fn handles_file_with_no_sequences() {
        let fasta_contents = ">seq1\n>seq2";
        let temp_file = create_temp_fasta_file(fasta_contents);
        let result = read_reference(temp_file.path());
        assert_eq!(result.len(), 2);
        assert_eq!(result["seq1"], Vec::<u8>::new());
        assert_eq!(result["seq2"], Vec::<u8>::new());
    }
}
