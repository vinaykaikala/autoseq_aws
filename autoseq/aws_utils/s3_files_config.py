#this module defines all the files required for each step


files_for_each_step = {
    'qc': {
        'base_dir': '/nfs/PROBIO/INBOX',
        'files': []
    },

    'alignment':{
        'base_dir': '/nfs/PROBIO/autoseq-genome',
        'files': [{'name': 'bwa', 'type': 'dir'}]
    }
}


files_for_each_step_new = {
    'qc': {
        'base': {'path': '/nfs/PROBIO/INBOX', 'files': []},
        'sample': {'path': '/nfs/PROBIO/autoseq-output', 'files': [] }
    }
}