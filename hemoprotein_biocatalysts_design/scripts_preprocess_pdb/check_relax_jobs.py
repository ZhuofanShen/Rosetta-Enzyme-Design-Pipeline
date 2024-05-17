import os
import shutil

os.chdir('HEM')
for pdb in os.listdir():
    if os.path.isdir(pdb + '/relax1'):
        err_file = None
        if os.path.isfile(pdb + '/relax1/' + pdb + '.err'):
            err_file = pdb + '/relax1/' + pdb + '.err'
        elif os.path.isfile(pdb + '/relax1/' + pdb + '1.err'):
            err_file = pdb + '/relax1/' + pdb + '1.err'
        if err_file:
            if not os.path.isfile(pdb + '/relax1/' + pdb + '_clean_relaxed.fasc'):
                with open(err_file, 'r') as pf:
                    lines = pf.readlines()
                if len(lines) > 0:
                    if lines[-1].endswith('out-of-memory handler.\n') or lines[-1].startswith('/var/lib/slurm/slurmd/job'):
                        shutil.rmtree(pdb + '/relax1')
                        if os.path.isdir(pdb + '/relax2'):
                            shutil.rmtree(pdb + '/relax2')
                        if os.path.isfile(pdb + '/3000'):
                            os.remove(pdb + '/3000')
                            with open(pdb + '/4000', 'w') as _:
                                pass
                        elif os.path.isfile(pdb + '/4000'):
                            os.remove(pdb + '/4000')
                            with open(pdb + '/5000', 'w') as _:
                                pass
                        elif os.path.isfile(pdb + '/5000'):
                            os.remove(pdb + '/5000')
                            with open(pdb + '/6000', 'w') as _:
                                pass
                        elif os.path.isfile(pdb + '/6000'):
                            os.remove(pdb + '/6000')
                            with open(pdb + '/7000', 'w') as _:
                                pass
                        elif os.path.isfile(pdb + '/7000'):
                            os.remove(pdb + '/7000')
                            with open(pdb + '/8000', 'w') as _:
                                pass
                        else:
                            with open(pdb + '/3000', 'w') as _:
                                pass
                        if lines[-1].startswith('/var/lib/slurm/slurmd/job'):
                            with open('../bus_error', 'a') as pf:
                                pf.write(pdb + '\n')
                    elif lines[-1].endswith('DUE TO PREEMPTION ***\n') or lines[-1].endswith('DUE TO TIME LIMIT ***\n'):
                        shutil.rmtree(pdb + '/relax1')
                        if os.path.isdir(pdb + '/relax2'):
                            shutil.rmtree(pdb + '/relax2')
                        with open(pdb + '/rerun', 'w') as _:
                            pass
                        with open('stucked', 'a') as pf:
                            pf.write(pdb + '\n')
                else:
                    shutil.rmtree(pdb + '/relax1')
                    if os.path.isdir(pdb + '/relax2'):
                        shutil.rmtree(pdb + '/relax2')
                    with open(pdb + '/rerun', 'w') as _:
                        pass
                    with open('../stucked', 'a') as pf:
                        pf.write(pdb + '\n')
            # else:
            #     for decoy in os.listdir(pdb + '/relax1'):
            #         if decoy.endswith('in_progress'):
            #             with open(pdb + '/continue', 'w') as _:
            #                 pass
        else:
            shutil.rmtree(pdb + '/relax1')
            if os.path.isdir(pdb + '/relax2'):
                shutil.rmtree(pdb + '/relax2')
            with open(pdb + '/rerun', 'w') as _:
                pass
    elif not os.path.isfile(pdb + '/' + pdb + '_relax.pdb') and not os.path.isfile(pdb + '/rerun') \
            and not os.path.isfile(pdb + '/3000') and not os.path.isfile(pdb + '/4000') \
            and not os.path.isfile(pdb + '/5000') and not os.path.isfile(pdb + '/6000') \
            and not os.path.isfile(pdb + '/no_proximal_res') and not os.path.isfile(pdb + '/myoglobin'):
        with open('../omitted', 'a') as pf:
            pf.write(pdb + '\n')
