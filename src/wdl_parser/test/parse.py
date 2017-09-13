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
            Int[Int]? num_references = 10000
            call PartyTime {
            '''

    # def test_match_keyword_argument(self):
    #     # todo improve tests
    #     data = 'Int? num_references = 10000'
    #     _ = wdl_parser.parse.VariableDefinitions.keyword_argument().parseString(data)

    def test_variable_definition(self):
        parsed = wdl_parser.parse.variable_definitions().parseString(self.data)
        self.assertEqual(len(parsed), 4)
        self.assertEqual(parsed[0]['variable_type'], 'File')
        self.assertEqual(parsed[1]['variable_type'], 'File')
        self.assertEqual(parsed[2]['variable_type'], 'String')

        self.assertEqual(parsed[0]['variable_name'], 'fasta')
        self.assertEqual(parsed[1]['variable_name'], 'gtf')
        self.assertEqual(parsed[2]['variable_name'], 'reference_name')

        self.assertEqual(parsed[3]['default_value'], '10000')
        self.assertEqual(parsed[3]['variable_name'], 'num_references')
        self.assertEqual(parsed[3]['variable_type'], 'Int[Int]?')


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
        parsed = wdl_parser.parse.call_input_assignment().parseString(self.data)['call_inputs']
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
        parsed = wdl_parser.parse.workflow_calls().parseString(self.data)['calls']

        self.assertEqual(len(parsed), 2)
        self.assertEqual(parsed[0]['task_name'], 'CellRangerMkref')
        self.assertEqual(len(parsed[0]['call_inputs']), 3)
        self.assertEqual(parsed[0]['call_inputs'][2]['variable_name'], 'reference_name')


class TestOutput(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = '''
            output {
              File reference_bundle = CellRangerMkref.reference_bundle
              File reference_bundle2 = CellRangerMkref.reference_bundle2
            }'''

        cls.data2 = '''
        output {
          Array[File] output_subset_fastqs = glob("./*_subset.fastq.gz")
        }'''

    def test_output(self):
        parsed = wdl_parser.parse.outputs().parseString(self.data)
        self.assertEqual(len(parsed['outputs']), 2)
        self.assertEqual(parsed['outputs'][1]['variable_type'], 'File')
        self.assertEqual(parsed['outputs'][1]['variable_name'], 'reference_bundle2')
        self.assertEqual(parsed['outputs'][1]['variable_value'],
                         'CellRangerMkref.reference_bundle2')

    def test_output2(self):
        parsed = wdl_parser.parse.outputs().parseString(self.data2)
        self.assertEqual(len(parsed['outputs']), 1)
        self.assertEqual(parsed['outputs'][0]['variable_type'], 'Array[File]')
        self.assertEqual(parsed['outputs'][0]['variable_name'], 'output_subset_fastqs')
        self.assertEqual(parsed['outputs'][0]['variable_value'], 'glob("./*_subset.fastq.gz")')



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
        # todo improve tests
        _ = wdl_parser.parse.workflow().parseString(self.data)


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

        cls.data2 = '''
        task SubsetFastqFromIndices {
          File indices_json
          File input_fastq_r1
          File input_fastq_r2
          File input_fastq_i1
        
          command <<<
            python3 <<CODE
        
            import os
            import json
            import scsequtil.fastq as fq
            import scsequtil.reader as rd
            import gzip
        
            # get indices
            with open('${indices_json}', 'r') as f:
              indices = set(json.load(f))
        
            # set fastq inputs
            fastqs = ['${input_fastq_r1}', '${input_fastq_r2}', '${input_fastq_i1}']
            readers = [fq.Reader(f) for f in fastqs]
        
            # define filenames
            filenames_nopath = [os.path.split(f)[1] for f in fastqs]
            output_filenames = [f.partition('.fastq')[0] + '_subset.fastq.gz' for f in filenames_nopath]
        
            # open some files
            output_fileobjs = [gzip.open(f, 'wt') for f in output_filenames]
        
            # write to file
            try:
                for records in rd.zip_readers(*readers, indices=indices):
                    for record, fout in zip(records, output_fileobjs):
                        fout.write(str(record))
            finally:
                for f in output_fileobjs:
                    f.close()
        
            CODE
          >>>
          runtime {
            docker: "ambrosejcarr/python3-scientific:0.2.0"
            memory: "2 GB"
            disks: "local-disk 220 HDD"
          }
          output {
            Array[File] output_subset_fastqs = glob("./*_subset.fastq.gz")
          }
        }'''

    def test_docker(self):
        data = 'docker: "humancellatlas/star_dev:v1"'
        parsed = wdl_parser.parse.docker().parseString(data)
        self.assertEqual(parsed['docker'], 'humancellatlas/star_dev:v1')

    def test_memory(self):
        data = 'memory: "8 GB"'
        parsed = wdl_parser.parse.memory().parseString(data)
        self.assertEqual(parsed['memory'], '8 GB')

    def test_disks(self):
        data = 'disks: "local-disk 220 HDD"'
        parsed = wdl_parser.parse.disks().parseString(data)
        self.assertEqual(parsed['disks']['disk_location'], 'local-disk')
        self.assertEqual(parsed['disks']['size'], '220')
        self.assertEqual(parsed['disks']['disk_type'], 'HDD')

    # @unittest.skip('not working yet')
    def test_parse_task(self):
        parsed = wdl_parser.parse.task().ignore(wdl_parser.parse.wdl_comment()
                                                ).parseString(self.data)
        print(parsed['task_name'])
        print(parsed['variable_definitions'])
        print(parsed['outputs'])

    # @unittest.skip('not working yet')
    def test_parse_task2(self):
        parsed = wdl_parser.parse.task().ignore(wdl_parser.parse.wdl_comment()
                                                ).parseString(self.data2)
        print(parsed)
        print(parsed['task_name'])
        print(parsed['variable_definitions'])
        print(parsed['outputs'])


if __name__ == "__main__":
    nose2.main()
