import unittest
import nose2
import wdl_parser.parse


class TestVariableDefinition(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        final line contains a badly spaced and formatted call that is technically
        correct and would validate.

        :return:
        """
        cls.data = '''
            File fasta
            File gtf
            String reference_name
            call PartyTime {
            '''

    def test_variable_definition(self):
        parsed = wdl_parser.parse._variable_definitions.parseString(self.data)
        self.assertEqual(len(parsed), 3)
        self.assertEqual(parsed[0]['type'], 'File')
        self.assertEqual(parsed[1]['type'], 'File')
        self.assertEqual(parsed[2]['type'], 'String')

        self.assertEqual(parsed[0]['variable_name'], 'fasta')
        self.assertEqual(parsed[1]['variable_name'], 'gtf')
        self.assertEqual(parsed[2]['variable_name'], 'reference_name')


class TestCallAssignedInput(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = '''
            input:
              fasta = fasta.weirdness,
              gtf = gtf,
              reference_name = i_am_different
            }'''

    def test_call_assigned_input(self):
        parsed = wdl_parser.parse._call_assigned_input.parseString(self.data)['inputs']
        self.assertEqual(len(parsed), 3)
        self.assertEqual(parsed[0]['value'], 'fasta.weirdness')
        self.assertEqual(parsed[1]['value'], 'gtf')
        self.assertEqual(parsed[2]['value'], 'i_am_different')

        self.assertEqual(parsed[0]['variable_name'], 'fasta')
        self.assertEqual(parsed[1]['variable_name'], 'gtf')
        self.assertEqual(parsed[2]['variable_name'], 'reference_name')


class TestCall(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = '''
            call CellRangerMkref {
              input:
                fasta = fasta,
                gtf = gtf,
                reference_name = reference_name
            }
            
            call CellRangerRandomThing {
              input:
                foo = bar
            }'''

    def test_call(self):
        parsed = wdl_parser.parse._calls.parseString(self.data)['calls']

        self.assertEqual(len(parsed), 2)
        self.assertEqual(parsed[0]['task_name'], 'CellRangerMkref')
        self.assertEqual(len(parsed[0]['inputs']), 3)
        self.assertEqual(parsed[0]['inputs'][2]['variable_name'], 'reference_name')


class TestOutput(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = '''
            output {
              File reference_bundle = CellRangerMkref.reference_bundle
              File reference_bundle2 = CellRangerMkref.reference_bundle2
            }'''

    def test_output(self):
        parsed = wdl_parser.parse._outputs.parseString(self.data)
        self.assertEqual(len(parsed['outputs']), 2)
        self.assertEqual(parsed['outputs'][1]['variable_type'], 'File')
        self.assertEqual(parsed['outputs'][1]['variable_name'], 'reference_bundle2')
        self.assertEqual(parsed['outputs'][1]['variable_value'],
                         'CellRangerMkref.reference_bundle2')


class TestWorkflow(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = '''
          workflow Generate10xDemoData {
            File fastq_r1
            File fastq_r2
            File fastq_i1
            File star_genome
            File gtf
            Int chromosome
        
            call StarAlignSubset {
              input:
                genomic_fastq = fastq_r2,
                gtf = gtf,
                star_genome = star_genome
            }
        
            call ExtractIndicesSpecificChromosomeAlignments {
              input:
                bam_file = StarAlignSubset.output_bam,
                chromosome = chromosome
            }
        
            call SubsetFastqFromIndices {
              input:
                indices_json = ExtractIndicesSpecificChromosomeAlignments.output_indices_json,
                input_fastq_r1 = fastq_r1,
                input_fastq_r2 = fastq_r2,
                input_fastq_i1 = fastq_i1
            }
        
            output {
              File subset_fastq_r1 = SubsetFastqFromIndices.output_subset_fastqs[0]
              File subset_fastq_r2 = SubsetFastqFromIndices.output_subset_fastqs[1]
              File subset_fastq_i1 = SubsetFastqFromIndices.output_subset_fastqs[2]
            }
          }'''

    def test_workflow(self):
        parsed = wdl_parser.parse._workflow.parseString(self.data)
        print(parsed)


class TestTask(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = '''
        task StarAlignSubset {
          File genomic_fastq
          File gtf
          File star_genome
          Int? subset_size = 1000000
        
          # note that STAR runThreadN must always equal 1 (default) or the order of the .bam file will be
          # disordered relative to the input data, breaking the expectation of same-ordering necessary
          # for downstream methods
          command {
            # unpack genome reference
            tar -zxvf ${star_genome}
        
            # truncate the input file
            zcat ${genomic_fastq} | head -n ${subset_size} > reads.fastq
        
            # align reads
            STAR  --readFilesIn reads.fastq \
              --genomeDir ./star \
              --outSAMstrandField intronMotif \
              --genomeLoad NoSharedMemory \
              --sjdbGTFfile ${gtf} \
              --outSAMtype BAM Unsorted  \
              --outSAMunmapped Within \
              --runThreadN 1
          }
          output {
            File output_bam = "Aligned.out.bam"
          }
          runtime {
            docker: "humancellatlas/star_dev:v1"
            memory: "8 GB"  # not used locally
            disks: "local-disk 220 HDD"  # 80G fastq, 32G reference bundle, 80G bam, 15% overflow (30G)
          }
        }'''

    def test_parse_task(self):
        print(self.data[124:135])
        parsed = wdl_parser.parse._task.parseString(self.data)


if __name__ == "__main__":
    nose2.main()
