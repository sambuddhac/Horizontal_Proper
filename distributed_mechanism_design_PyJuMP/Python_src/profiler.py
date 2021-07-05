from time import time

class Profiler(object):

   def __init__(self):
      self.start_time = None
      self.end_time = None
      self.duration = None
      self.last_time = None

   def start(self):
      self.start = time()

   def stop(self):
      self.end = time()

      if self.start is not None:
          self.duration = self.end - self.start

   def get_duration(self):
       return self.duration

   def get_elapsed(self):
       if self.start is not None:
           return time() - self.start
       return "profiler not started"

   def get_interval(self):
       if self.start is not None and self.last_time is None:
           self.last_time = time()
           return self.last_time - self.start
       elif self.start is not None:
           delta = time() - self.last_time
           self.last_time = time()
           return delta
       return "profiler not started"
