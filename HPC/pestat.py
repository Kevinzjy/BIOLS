#! /usr/bin/python3
# -*- encoding:utf-8 -*-
# modified from pydf [https://github.com/garabik/pydf]
# Usage:
#   pestat all/high/GPUfat/dev/middle
import sys
import unicodedata
from collections import defaultdict

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

try:
    from commands import  getstatusoutput
except ImportError:
    from subprocess import getstatusoutput


# Home of PBS directory
#PBS_HOME='/usr/local/torque-4.2.9/bin/'
PBS_HOME='/usr/local/bin/'

# Width of each row
FORMAT = [
    ('node', 10, "l"), ('state', 7, "r"), ('load', 6, "r"),
    ('used', 7, "r"), ('mem', 7, "r"), ('avail', 7, "r"), ('ncpus', 5, "r"), ('perc', 7, "r"),
    ('bar', 30, "l"), ('queue', 11, "l"), ('job_list', 11, 'l'),
]

# Exclude queue in blacklist
QUEUE_BLACKLIST = ['manager']

#some default definitions
COLORS = {
    'none'       :    "",
    'default'    :    "\033[0m",
    'bold'       :    "\033[1m",
    'underline'  :    "\033[4m",
    'blink'      :    "\033[5m",
    'reverse'    :    "\033[7m",
    'concealed'  :    "\033[8m",

    'black'      :    "\033[30m",
    'red'        :    "\033[31m",
    'green'      :    "\033[32m",
    'yellow'     :    "\033[33m",
    'blue'       :    "\033[34m",
    'magenta'    :    "\033[35m",
    'cyan'       :    "\033[36m",
    'white'      :    "\033[37m",

    'on_black'   :    "\033[40m",
    'on_red'     :    "\033[41m",
    'on_green'   :    "\033[42m",
    'on_yellow'  :    "\033[43m",
    'on_blue'    :    "\033[44m",
    'on_magenta' :    "\033[45m",
    'on_cyan'    :    "\033[46m",
    'on_white'   :    "\033[47m",

    'beep'       :    "\007"
}


BAR = '#'
FILL = '.'
class Bar(object):
    "Progress bar"
    def __init__(self, percentage=0, width=30, header=False):
        self.percentage = percentage
        self.width = width
        self.header = header

    def __len__(self):
        return self.width

    def __str__(self):
        return self.format('l')

    def format(self, pos):
        if self.header:
            return ' '*self.width
        size = min(int(round(self.percentage*(self.width-2))), self.width-2)
        return '[' + manglestring(size * BAR, self.width - 2, pos, FILL) + ']'


STATE_COLOR = {
    'free': "green",
    'busy': "cyan",
    'excl': "red",
    'down': "default",
    'offl': "on_red",
    'UNKN': "default",
}

STATE = {
    'job-exclusive': 'excl',
    'job-exclusive,busy': 'busy',
    'busy': 'busy',
    'free': 'free',
    'offline': 'offl',
    'offline,job-exclusive': 'offl',
    'offline,job-exclusive,busy': 'offl',
    'down': 'down',
    'down,busy': 'down',
    'down,offline': 'offl',
    'down,job-exclusive': 'down',
    'down,offline,job-exclusive': 'offl',
    'down,offline,busy': 'offl',
    'down,offline,job-exclusive, busy': 'offl',
    'UNKN': 'UNKN',
}

def get_stat(job_info, target_queue):
    """Parse pbsnodes output"""
    nodes_stat = []
    status, output = getstatusoutput(PBS_HOME + "pbsnodes -aox")
    root = ET.fromstring(output)
    for node in root:
        name = node.find('name').text
        queue = node.find('properties').text
        if queue in QUEUE_BLACKLIST:
            continue

        if target_queue != "all" and queue != target_queue:
            continue

        state = STATE[node.find('state').text]

        ncpus = int(node.find('np').text)

        # optional field
        status = node.find('status')
        if status is None:
            mem = None
            aval = None
            used = None
            loadave = None
            msg = None
        else:
            status_d = dict((i.split('=')[0], i.split('=')[1]) for i in status.text.split(','))

            mem = float(status_d['physmem'][:-2]) / 1024 / 1024
            tmp_mem = float(status_d['totmem'][:-2]) / 1024 / 1024
            tmp_aval = float(status_d['availmem'][:-2]) / 1024 / 1024
            used = max(tmp_mem - tmp_aval, 0)
            aval = max(mem - used, 0)
            loadave = float(status_d['loadave'])
            
            msg = status_d['message'] if 'message' in status_d else None

        # load = float(status_d['loadave'])
        jobs = node.find('jobs')

        tmp_jobs = {}
        core = 0
        if jobs is not None:
            node_jobs = jobs.text.split(',')
            for i in node_jobs:
                job_cpu, job_id = i.split('/')
                job_cores = job_cpu.split('-')
                if len(job_cores) > 1:
                    tmp_jobs[job_id] = tmp_jobs.setdefault(job_id, 0) + int(job_cores[1]) - int(job_cores[0]) + 1
                else:
                    tmp_jobs[job_id] = tmp_jobs.setdefault(job_id, 0) + 1
            core += sum([tmp_jobs[i] for i in tmp_jobs])

        if core == 0 and state == 'busy':
            state = 'UNKN'

        # User on node
        user_on_node = []
        if tmp_jobs:
            for i in tmp_jobs:
                if i in job_info:
                    user_on_node.append('{0}/{1} {2}'.format(tmp_jobs[i], i.split('.')[0], job_info[i]))
                else:
                    continue
                    # user_on_node.append('{0}/{1}'.format(tmp_jobs[i], i.split('.')[0]))
        user_str = ', '.join(user_on_node)

        node_info = [
            (COLORS['default'], name),
            (COLORS[STATE_COLOR[state]], state),
            (COLORS[STATE_COLOR[state]], loadave),
            (COLORS[STATE_COLOR[state]], '' if used is None else '{0:.1f}'.format(used)),
            (COLORS[STATE_COLOR[state]], '' if mem is None else '{0:.1f}'.format(mem)),
            (COLORS[STATE_COLOR[state]], '' if aval is None else '{0:.1f}'.format(aval)),
            (COLORS['default'], ncpus),
            (COLORS[STATE_COLOR[state]], '{0:2d}'.format(core)),
            (COLORS[STATE_COLOR[state]], Bar(percentage=float(core)/ncpus)),
            (COLORS['default'], queue),
            (COLORS['default'], user_str if msg is None else msg),
        ]
        nodes_stat.append(node_info)
    return nodes_stat


def manglestring(s, l, pos, fillchar=' '):
    "cut string to fit exactly into l chars"
    if pos == "r":
        ns = str.rjust(s, l, fillchar)
    elif pos == "l":
        ns = str.ljust(s, l, fillchar)
    elif pos == "c":
        ns = str.center(s, l, fillchar)
    else:
        raise ValueError('Error in manglestring')
    if len(ns) > l:
        ns = ns[:int(l/2)] + "~" + ns[-int(l/2)+1:]
    return ns


def out(s):
    try:
        sys.stdout.write(s)
    except UnicodeEncodeError:
        sys.stdout.write(s.encode('ascii', 'ignore').decode())


HEADER = {
    'node': "Node",
    'state': "State",
    'load': "Load%",
    'used': "Used",
    'mem': "Mem",
    'avail': "Avail",
    'ncpus': 'CPU',
    'perc': "Proc",
    'bar': str(Bar(header=True)),
    'queue': "Queue",
    'job_list': "Job_List",
}
def display_table(table):
    header_colour = "yellow"
    for i, j in enumerate(FORMAT):
        width, pos = FORMAT[i][1:3]
        out(COLORS[header_colour])
        out(manglestring(HEADER[str(j[0])], width, pos))
        out(COLORS[header_colour])
        out(' ')
    out('\n')

    for row in table:
        for i, j in enumerate(row[:-1]):
            width, pos = FORMAT[i][1:3]
            out(j[0])
            out(manglestring(str(j[1]), width, pos))
            out(COLORS['none'])
            out(' ')
        out(COLORS['default'])
        out(row[-1][1])
        out(COLORS['default'])
        out('\n')


def get_jobinfo():
    """Get username for each job"""
    job_info = {}
    status, output = getstatusoutput(PBS_HOME + 'qstat -r')
    if status != 0:
        return None

    if output:
        for line in output.splitlines()[5:]:
            content = line.split()
            job_id, username = content[0], content[1]
            job_info[job_id] = username
    return job_info


def main():
    """Main func"""
    queue = "all"
    if len(sys.argv) > 1:
        queue = sys.argv[1]

    job_info = get_jobinfo()
    if job_info is None:
        sys.exit('Unable to communicate with PBS server host')
    nodes_info = get_stat(job_info, queue)
    display_table(nodes_info)


if __name__ == "__main__":
    main()
