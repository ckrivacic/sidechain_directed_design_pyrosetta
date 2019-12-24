import sys, os, subprocess, re
'''
Input mover type as first argument
'''
'''
Please put path to pdbredo as $PDBREDO in your bashrc
Usage:
    test_movers.py <input_df> <output_dir> <mover> [br]
'''
task_len = 500

def file_len(fname):
    with open(fname) as f:
        for i,l in enumerate(f):
            pass
    return i + 1

def submit(alignments, **params):
    from klab import cluster, process
    import json
    csv_path = sys.argv[1]
    #num_tasks = (file_len(alignments) // task_len) + 1
    output_directory = 'logs'
    #error_directory = 'errors'
    mover = str(sys.argv[3])
    outdir = str(sys.argv[2])
    if len(sys.argv) > 4:
        backrub = str(sys.argv[4])
    else:
        backrub = ''
    if backrub=='br':
        #num_tasks = file_len(csv_path) * task_len
        num_tasks = 100 * task_len
    else:
        num_tasks = 14 * task_len

    max_runtime = params.get('max_runtime','12:00:00')
    max_memory = params.get('max_memory','4G')

    script_path = os.path.expanduser('~/intelligent_design/sidechain_directed_design_pyrosetta/inverse_rotamers/test_movers.py')

    qsub_command = 'qsub', '-h', '-cwd',
    qsub_command += '-o', output_directory,
    #qsub_command += '-e', error_directory,
    qsub_command += '-j', 'y',
    qsub_command += '-t', '1-{0}'.format(num_tasks),
    qsub_command += '-l', 'h_rt={0}'.format(max_runtime),
    qsub_command += '-l', 'mem_free={0}'.format(max_memory),
    qsub_command += script_path,
    qsub_command += csv_path,
    qsub_command += outdir,
    qsub_command += mover,
    qsub_command += backrub,
    print(qsub_command)

    status = process.check_output(qsub_command, stderr=subprocess.STDOUT).decode('utf-8')

    status_pattern = re.compile(r'Your job-array (\d+).[0-9:-]+ \(".*"\) has been submitted')
    status_match = status_pattern.match(status)

    if not status_match:
        print(status)
        sys.exit()

    # Figure out the job id, then make a params file for it.
    job_id = status_match.group(1)

    qrls_command = 'qrls', job_id
    process.check_output(qrls_command)
    print(status)

if __name__ == '__main__':
    submit('alignment_results.m8')