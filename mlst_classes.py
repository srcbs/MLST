# class
class Moab:
   '''Submits a list of calls to the scheduler using msub/xmsub. Job dependencies are controlled using depend (logical), depend_type ('one2one', 'expand', 'conc', 'complex', 'all'), depend_val (list of integers) and ids (list of jobids):
         
         'one2one' makes job 1 dependent on id 1, job 2 on id 2 ...
         'expand' makes n first jobs dependent on id 1, next n jobs dependent on id 2 ...
         'conc' makes job 1 dependent on n first ids, job 2 dependent on next n ids ...
         'complex' takes several integers in a list as input and makes the jobs dependent on the number of ids given. Eg.[2,1] means that the first job will be dependent on the first 2 ids and the second job will be dependent on the 3 id ...
      
      Jobs are submitted as hold by default and should be released using Moab.release().
   '''
   
   def __init__(self, calls, logfile=None, runname='run_test', queue='cbs', cpu='nodes=1:ppn=1,mem=2gb,walltime=43200', depend=False, depend_type='one2one', depend_val=[], hold=True, depend_ids=[], partition='cge-cluster', env=None, host=None):
      '''Constructor for Moab class'''
      self.calls = calls
      self.runname = runname
      self.logfile = logfile
      self.queue = queue
      self.cpu = cpu
      self.depend = depend
      self.depend_type = depend_type
      self.depend_val = depend_val
      self.hold = hold
      self.depend_ids = depend_ids
      self.env = env
      self.partition = partition
      self.host = host
      
      # put jobs in queue upon construction
      self.dispatch()
   
   def __repr__(self):
      '''Return string of attributes'''
      msg = 'Moab(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)' % ("["+", ".join(self.calls)+"]", self.logfile, self.runname, self.queue,  self.cpu, str(self.depend), self.depend_type, "["+", ".join(map(str,self.depend_val))+"]", str(self.hold), "["+", ".join(self.depend_ids)+"]", self.env, self.partition, self.host)
      return msg
   
   def get_logger(self):
      '''Return logger object'''
      
      if self.logfile:
         import logging
         
         # if already a logger do nothing and return
         # else create logger with given logfile
         if isinstance(self.logfile, logging.Logger):
            return self.logfile
         else:
            logger = logging.getLogger('moab_submissions')
            hdlr = logging.FileHandler('%s' % self.logfile)
            formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
            hdlr.setFormatter(formatter)
            logger.addHandler(hdlr) 
            logger.setLevel(logging.INFO)
            return logger
      else:
         return None
   
   def create_dependencies(self):
      '''Create job dependencies
         
         'one2one' makes job 1 dependent on id 1, job 2 on id 2 ...
         'expand' makes n first jobs dependent on id 1, next n jobs dependent on id 2 ...
         'conc' makes job 1 dependent on n first ids, job 2 dependent on next n ids ...
         'complex' takes a n=list as input and makes the jobs dependent on the number of ids given in n. Eg.[2,1] means that the first job will be dependent on the first 2 ids and the second job will be dependent on the 3 id ...
      '''
      
      if not self.depend:
         depends = None
      else:
            if self.depend_type == 'one2one':
               depends = self.depend_ids
            
            elif self.depend_type == 'expand':
               n = int(self.depend_val[0])
               depends = []
               depend_ids = self.depend_ids     # do not want to remove from self.depend_ids because it will remove the ids from input class instance (eg. other calls)
               c = 0
               for j in range(len(self.calls)):
                  c = c + 1
                  if c < int(n):
                     depends.append(depend_ids[0])
                  if c == int(n):
                     depends.append(depend_ids[0])
                     depend_ids = depend_ids[1:]
                     c = 0
            
            elif self.depend_type == 'conc':
               n = int(self.depend_val[0])
               depends = []
               depend_ids = self.depend_ids     # do not want to remove from self.depend_ids because it will remove the ids from input class                
               for j in range(len(self.calls)):
                  s = ':'.join(depend_ids[:int(n)])
                  depends.append(s)
                  depend_ids = depend_ids[int(n):]
            
            elif self.depend_type == 'complex':
               old_index = 0
               depends = []
               for index in self.depend_val:
                  s = ':'.join(self.depend_ids[old_index:(index+old_index)])
                  depends.append(s)
                  old_index=index
            
            elif self.depend_type == 'all':
               depends = []
               for j in range(len(self.calls)):
                  curr = ':'.join(self.depend_ids)
                  depends.append(curr)
            
            else:
               raise AttributeError('depend_type not recognized: %s' % self.depend_type)
      return depends
   
   def ssh_submit(self, host, cmd):
      '''Will ssh to host, submit commands and return job ids'''
      
      import paramiko, base64
      import os
      
      client = paramiko.SSHClient()
      client.load_system_host_keys()
      client.set_missing_host_key_policy(paramiko.WarningPolicy())
      client.connect(host)
      
      stdin, stdout, stderr = client.exec_command(cmd)
      for line in stdout:
         id = line.rstrip('\n')
      err_msg = ''.join(stderr)
      
      client.close()      
      return id, err_msg

   def submit_xmsub(self, depends, logger):
      '''Submits jobs using xmsub'''
      
      import re
      import subprocess
      import time
      import os
      import mlst_modules
      
      home = os.getcwd()
      paths = mlst_modules.setSystem()
      
      ids = []
      for i in range(len(self.calls)):
         call = self.calls[i]
         stdout = '%s/log/%s%i.o' % (home, self.runname, i)
         stderr = '%s/log/%s%i.e' % (home, self.runname, i)
         
         # catch stdouts if call includes 'program infile > outfile', needs to be directed as -O instead of >
         pattern = re.compile(r'(^.+)>\s(.+)$')
         match = pattern.search(call)
         if match:
            call = match.group(1)
            stdout = '%s/%s' % (home, match.group(2))
         
         # create xmsub commands
         cmd = paths['mlst_home'] + 'xmsub'
         
         # toggle if job should be on hold or env variable should be added
         if self.hold: cmd = '%s -h ' % cmd
         if self.env: cmd = cmd + ' -v %s' % self.env
         
         if not self.depend:
            xmsub = cmd+' -d %s -l %s,partition=%s -O %s -E %s -r y -q %s -N %s -t %s' % (home, self.cpu, self.partition, stdout, stderr, self.queue, self.runname, call)
         else:
            xmsub = cmd+' -d %s -l %s,depend=%s,partition=%s -O %s -E %s -r y -q %s -N %s -t %s' % (home, self.cpu, depends[i], self.partition, stdout, stderr, self.queue, self.runname, call)
         
         time.sleep(1)
         if logger: logger.info(xmsub)
         
         # submit on different host if that is given
         if self.host:
            try:
               (id, stderr) = self.ssh_submit(self.host, xmsub)
            except:
               print stderr
               print 'Job error, waiting 1m'
               time.sleep(60)
               (id, stderr) = self.ssh_submit(self.host, xmsub)
            ids.append(id)
         else:
            try:
               id = subprocess.check_output(xmsub, shell=True)
            except:
               print 'Job error, waiting 1m'
               time.sleep(60)
               id = subprocess.check_output(xmsub, shell=True)
            ids.append(id.split('\n')[1])
      return ids
   
   def submit_wrapcmd(self, depends, logger):
      '''Take input as command and submit to msub (this way pipes and redirects can be done)'''
      
      import subprocess
      import random
      import string
      import time
      import os
      
      home = os.getcwd()
            
      ids = []
      for i in range(len(self.calls)):
         call = self.calls[i]
         N = 10
         rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(N))
         filename = '%s/pbsjob.tmp%s' % (home, rand)
         
         stdout = '%s/log/%s%i.o' % (home, self.runname, i)
         stderr = '%s/log/%s%i.e' % (home, self.runname, i)
         
         # write pbsjob file
         fh = open(filename, 'w')
         fh.write('#!/bin/sh\n\n')
         #fh.write('newgrp cge\n')      # temporary fix to group problems on cge-s2 cluster
         fh.write('%s\n' % call)
         fh.close()
         
         # create msub command
         cmd = 'msub'
         # toggle if on hold or env variable
         if self.hold: cmd = '%s -h' % cmd
         if self.env: cmd = cmd + ' -v %s' % self.env
         
         if not self.depend:
            msub = '%s -d %s -l %s,partition=%s -o %s -e %s -q %s -r y -N %s %s' % (cmd, home, self.cpu, self.partition, stdout, stderr, self.queue, self.runname, filename)
         else:
            msub = '%s -d %s -l %s,depend=%s,partition=%s -o %s -e %s -q %s -r y -N %s %s' % (cmd, home, self.cpu, depends[i], self.partition, stdout, stderr, self.queue, self.runname, filename)
         
         time.sleep(1)
         if logger: logger.info(msub)

         # submit on different host if that is given
         if self.host:
            try:
               (id, stderr) = self.ssh_submit(self.host, msub)
            except:
               print stderr
               print 'Job error, waiting 1m'
               time.sleep(60)
               (id, stderr) = self.ssh_submit(self.host, msub)
            ids.append(id)
         else:
            try:
               id = subprocess.check_output(msub, shell=True)
            except:
               print 'Job error, waiting 1m'
               time.sleep(60)
               id = subprocess.check_output(msub, shell=True)
            ids.append(id.split('\n')[1])
         
         # remove pbsjob file
         #rm_files([filename])
         #print ids
      return ids
   
   def dispatch(self):
      '''Submit job to queue, if depend=True dependent jobids should be given as ids. 
         Dependtype can be 'one2one', 'expand' 'conc', 'complex'.
         
         'one2one' makes job 1 dependent on id 1, job 2 on id 2 ...
         'expand' makes n first jobs dependent on id 1, next n jobs dependent on id 2 ...
         'conc' makes job 1 dependent on n first ids, job 2 dependent on next n ids ...
         'complex' takes a n=list as input and makes the jobs dependent on the number of ids given in n. Eg.[2,1] means that the first job will be dependent on the first 2 ids and the second job will be dependent on the 3 id ...'''
      
      # start logger
      logger = self.get_logger()
      
      # create job dependencies
      depends = self.create_dependencies()
      
      if self.calls[0].find('|') > -1 or self.calls[0].find('<') > -1 or self.calls[0].find('>>') > -1 or self.calls[0].find(';') > -1:
         # perform wrapcmd if calls includes pipes / left-redirects
         self.ids = self.submit_wrapcmd(depends, logger)
      else:
         # perform xmsub if calls does not include pipes (can have right-redirects)
         self.ids = self.submit_xmsub(depends, logger)
      
   
   def release(self, host):
      '''Release submitted jobs from hold'''
      
      import subprocess
      
      if host:
         import paramiko, base64
         import os
         
         client = paramiko.SSHClient()
         client.load_system_host_keys()
         client.set_missing_host_key_policy(paramiko.WarningPolicy())
         client.connect(host)
         
         while len(self.ids) > 0:
            if len(self.ids) > 199:
               cmd = 'mjobctl -u user \"%s\"' % (' ').join(self.ids[:199])
               del self.ids[:199]
               stdin, stdout, stderr = client.exec_command(cmd)
            else:
               cmd = 'mjobctl -u user \"%s\"' % (' ').join(self.ids)
               stdin, stdout, stderr = client.exec_command(cmd)
               break
      else:
         while len(self.ids) > 0:
            if len(self.ids) > 199:
               cmd = 'mjobctl -u user \"%s\"' % (' ').join(self.ids[:199])
               del self.ids[:199]
               out = subprocess.check_output(cmd, shell=True)
            else:
               cmd = 'mjobctl -u user \"%s\"' % (' ').join(self.ids)
               out = subprocess.check_output(cmd, shell=True)
               break


class Semaphore:
   '''Wait for files to be created, times are in seconds'''
   
   def __init__(self, semaphore_ids, home, file_prefix, queue, check_interval, max_time, partition, host=None):
      '''Constructor for Semaphore class'''
      self.semaphore_ids = semaphore_ids
      self.home = home
      self.file_prefix = file_prefix
      self.queue = queue
      self.check_interval = check_interval
      self.max_time = max_time
      self.host = host
   
   def ssh_submit(self, host, cmd):
      '''Will ssh to host, submit commands and return job ids'''
      
      import paramiko, base64
      import os
      
      client = paramiko.SSHClient()
      client.load_system_host_keys()
      client.set_missing_host_key_policy(paramiko.WarningPolicy())
      client.connect(host)
      
      stdin, stdout, stderr = client.exec_command(cmd)
      id = ''
      for line in stdout:
         id = line.rstrip('\n')
      err_msg = ''.join(stderr)
      
      client.close()      
      return id, err_msg
   
   def wait(self):
      '''Wait for files to be created'''
      
      from time import sleep
      import string
      import random
      import os
      import mlst_modules
      import subprocess
      
      paths = mlst_modules.setSystem()
      
      # add directory and set semaphore filename
      if not os.path.exists('semaphores/'):
         os.makedirs('semaphores/')
      
      rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(10))
      semaphore_file = 'semaphores/' + self.file_prefix + '.' + rand
      semaphore_file_err = 'log/' + self.file_prefix + '.' + rand + '.err'
      
      # create job 
      depends = ':'.join(self.semaphore_ids)
      xmsub = '%sxmsub -d %s -l ncpus=1,mem=10mb,walltime=180,depend=%s,partition=%s -O %s -q %s -N semaphores -E %s -r y -t echo done' % (paths['mlst_home'], self.home, depends, partition, semaphore_file, self.queue, semaphore_file_err)
      
      # submit job
      if self.host:
         dummy_id, stderr = self.ssh_submit(self.host, xmsub)
         if stderr:
            print stderr
      else:
         dummy_id = subprocess.check_output(xmsub, shell=True)
      
      # check for file to appear
      cnt = self.max_time
      while cnt > 0:
         if os.path.isfile(semaphore_file):
            break
         cnt -= self.check_interval
         sleep(self.check_interval)
      if cnt <= 0:
         raise SystemExit('%s did not finish in %is' % ())

class Semaphore_file:
   '''Wait for a file to appear'''
   
   def __init__(self, semaphore_file, check_interval, max_time):
      '''Constructor for Semaphore_file class'''
      self.semaphore_file = semaphore_file
      self.check_interval = check_interval
      self.max_time = max_time
   
   def wait(self):
      '''Wait for files to be created'''
      
      import os
      import subprocess
      from time import sleep
      
      # add directory and set semaphore filename
      if not os.path.exists('semaphores/'):
         os.makedirs('semaphores/')
      
      # check for file to appear
      cnt = self.max_time
      while cnt > 0:
         if os.path.isfile(self.semaphore_file):
            break
         cnt -= self.check_interval
         sleep(self.check_interval)
      if cnt <= 0:
         raise SystemExit('%s did not finish in %is' % ())
   
   
   