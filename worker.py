#!/usr/bin/env python3
import argparse
import requests
import sys
import subprocess
import datetime
import pytimeparse
import logging


parser = argparse.ArgumentParser(description='Schedules tasks')
parser.add_argument('-exp', '--exp_id', type=int, help='Experiment ID', nargs='?', default=None)
parser.add_argument('-exe', '--executable', type=str, help='Executable file', nargs='?', default=None)

parser.add_argument('-task', '--task_id', type=int, help='Task ID', nargs=None, default=None)
parser.add_argument('-calls', '--calls', type=int, help='Number of calls', nargs=None, default=None)
parser.add_argument('-duration', '--duration', type=float, help='Execution duration', nargs=None, default=None)
parser.add_argument('-subs', '--subregions', type=int, help='Number of subregions', nargs=None, default=None)
parser.add_argument('-st', '--status', type=str, help='Status of the task. D - done, S - suspended.', nargs=None, default=None)


def run_next_task(exp_id, executable):
    url = 'http://dakis.gimbutas.lt/api/exp/%d/next-task/' % exp_id
    resp = requests.get(url)
    task = resp.json()
    logging.info('Getting next task: exp_id=%d, status_code=%d, resp=%s' % (exp_id, resp.status_code, resp.json()))
    if resp.status_code == 200 and resp.json():
        cmd = [
            './%s' % executable.strip('./'),
            '--gkls_cls=%d' % task['func_cls'],
            '--gkls_fid=%d' % task['func_id'],
            '--task_id=%d' % task['id'],
            '--callback=%s' % sys.argv[0].strip('./'),
        ]
        logging.info('Calling: cmd=%s' % ' '.join(cmd))
        try:
            subprocess.Popen(cmd, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
        except Exception as e:
            logging.error('Got error while calling:  %s' % str(e))
        # logging.info(': cmd=%s' % ' '.join(cmd))
        # print('Started', ' '.join(cmd))
        # stdin=None, stdout=None, stderr=None, close_fds=True)


def send_task_results(args):
    url = 'http://dakis.gimbutas.lt/api/tasks/%d/' % args.task_id
    data = {
        'calls': args.calls,
        'duration': args.duration,
        'subregions': args.subregions,
        'status': args.status,
    }
    logging.info('Sending: url=%s, data=%s' % (url, data))
    resp = requests.put(url, data)
    logging.info('Response from sending: status_code=%s, data=%s' % (resp.status_code, resp.json()))
    exp_id = resp.json()['experiment'].split('experiments')[-1].strip('/')
    return int(exp_id)


def main(args):
    logging.basicConfig(filename='worker.log', level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',)
    logging.info('Invoked with: exp_id=%s executable=%s task_id=%s calls=%s '
                 '-duration=%s --subregions=%s --status=%s' % (args.exp_id, args.executable,
                 args.task_id, args.calls, args.duration, args.subregions, args.status))
    if args.exp_id:
        if not args.executable:
            raise ValueError('Provide executable file (-exe=) for this experiment')
        run_next_task(args.exp_id, args.executable)
    elif args.task_id:
        exp_id = send_task_results(args)
        run_next_task(exp_id, args.executable)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
