import os
from sarge import capture_stdout


def _run(cmd):
    return capture_stdout(cmd).stdout.text


def _normalize_names(fn, fn_new):
    with open(fn) as fh:
        with open(fn_new, 'w+') as fh_new:
            for l in fh:
                l = l.strip().split(' ', 1)[0]
                fh_new.write('{}\n'.format(l))


def _handle_read(base_path, base_name, read, new_path, ext='.trim.fasta'):
    _run('gzip -k -d {}/{}_{}_001{}.gz'.format(base_path, base_name,
                                               read, ext))
    _normalize_names(
        '{}/{}_{}_001{}'.format(base_path, base_name, read, ext),
        '{}/{}_{}_001.sync{}'.format(new_path, base_name, read, ext))


def run_align(args):
    print 'Starting alignment'
    ext = args.ext
    for study in os.listdir(args.base_dir):
        if args.studies is not None and study not in args.studies:
            print 'SKIPPING {}'.format(study)
            continue
        print '\tIn study {}'.format(study)
        for date in os.listdir('/'.join([args.base_dir, study])):
            print '\t\t{}'.format(date)
            for sample in os.listdir('/'.join(
                    [args.base_dir, study, date, 'raw'])):
                base_name = sample.rsplit('.', 1)[0].rsplit('_', 2)[0]
                base_path = '/'.join([args.base_dir, study, date, 'raw'])
                new_path = '/'.join([args.base_dir, study, date, 'processed'])
                presto_path = '/'.join([args.base_dir, study, date, 'presto'])
                if os.path.isfile('{}/{}.log'.format(presto_path, base_name)):
                    print ('Skipping {} since it already '
                           'exists.').format(base_name)
                    continue

                try:
                    os.makedirs(new_path)
                    os.makedirs(presto_path)
                except:
                    pass

                print '\t\t\t{}'.format(base_name)
                _handle_read(base_path, base_name, 'R1', new_path, ext)
                if not os.path.isfile('{}/{}_{}_001{}.gz'.format(
                        base_path, base_name, 'R2', ext)):
                    print '\t\t\t\tCan\'t find R2.  Using R1 alone.'
                    continue
                else:
                    _handle_read(base_path, base_name, 'R2', new_path, ext)

                try:
                    os.remove('{}/{}_{}_001{}'.format(base_path,
                              base_name, 'R1', ext))
                    os.remove('{}/{}_{}_001'.format(base_path,
                              base_name, 'R2', ext))
                except:
                    pass
                print _run(('AssemblePairs.py align -1 {0}/{1}_R1_001.sync{3} '
                            '-2 {0}/{1}_R2_001.sync{3} --fasta '
                            '--rc tail --log {2}/{1}.log --outdir {2} --nproc '
                            '4').format(new_path, base_name, presto_path, ext))
