#!/usr/bin/env python3

from collections import deque, namedtuple
from os import stat

def highest_power_of_two(n):
  if n < 0:
    raise ValueError('positive values only')
  count = 0
  n = n >> 1
  while n:
    n = n >> 1
    count += 1
  return count



class DecimationBuffer(deque):
  def __init__(self, n):
    super().__init__([0] * n, n)

  def feed(self, data):
    self.extend(data)



class Decimator:
  def __init__(self, n):
    self._buflen = n
    self._max_bufs = (highest_power_of_two(n) + 1)
    self._bufs = [None] * self._max_bufs
    self.decimate(0)

  def __getitem__(self, idx):
    return self._bufs[idx]
  
  def __iter__(self):
    return ((pow, b) for pow, b in enumerate(self._bufs) if b is not None)

  def decimate(self, pow):
    if pow >= self._max_bufs:
      raise ValueError(f'maximum decimation power is lg(n)')
    if self._bufs[pow] is None:
      self._bufs[pow] = DecimationBuffer(self._buflen)

  def feed(self, data):
    for pow, buf in enumerate(self._bufs):
      if buf is not None:
        if pow == 0:
          buf.feed(data)
        else:
          skips = 0x01 << pow
          buf.feed(data[::skips])
  
class FFTQuery:
  Segment = namedtuple('Segment', ['buf', 'width', 'domain'])

  @staticmethod
  def _closest_idx(f, width, len):
    l = int(f // width)  # low index
    lf = l * width  # corresponding freq
    
    idx = l         # assume low index

    if idx < (len - 1):  # only in this case can we offer a higher segment
      diff = abs(f - lf)          # check the difference from low freq
      if diff > (width/2):        # if greater than half width
        idx += 1                  # add one to index
    
    return idx

  @staticmethod
  def _bookend_idxs(f, width, len):
    l = int(f // width)  # low index
    idx = l         # assume low index
    if idx > (len - 2):
      raise Exception('invalid bookend indices')
    return (idx, idx + 1)

  @staticmethod
  def _floor_idx(f, width):
    return int(f // width)

  def __init__(self, n, fs):
    self._base = fs/n
    self._seglen = n//2
    self._segments = [None] * (highest_power_of_two(n) + 1)
    self.track(0)

  def _seg_by_index(self, idx):
    # since segments maintain computed index domains we can (relatively) quickly find the proper data
    pow = len(self._segments)
    seg = None
    while seg is None:
      pow -= 1
      if pow < 0:
        break
      s = self._segments[pow]
      if s is None:
        continue
      # check whether segment s contains idx:
      if s.domain['low'] <= idx and s.domain['high'] >= idx:
        seg = s
    if seg is None:
      # in this case no segment contianed the desired index
      raise IndexError()
    return seg

  def _internal__setitem__(self, idx, val):
    seg = self._seg_by_index(idx)
    b, _, domain = seg
    segindex = idx - domain['low'] + domain['covered']
    b[segindex] = val
  
  def __getitem__(self, idx):
    seg = self._seg_by_index(idx)  
    b, width, domain = seg  # get desired entry from this segment
    segindex = idx - domain['low'] + domain['covered']  # convert desired index into segment index
    return (segindex * width, b[segindex])
  
  def _recompute_domains(self):
    pow = len(self._segments)

    start_idx = 0
    fmax = 0      # maximum frequency covered so far

    while True:
      pow -= 1
      if pow < 0: # done when pow goes negative
        break

      seg = self._segments[pow]
      if seg is None:
        continue

      # how many indices can this segment contribute?
      _, width, domain = seg

      num_covered = int(fmax // width) + 1  # fmax has already been accounted for, and fmax // width is are the indices already covered
                                            # note: we add +1 to num_covered to eliminate the dc component of all segments
      num_contibuting = int(self._seglen - num_covered) # this segment will provide any indices that haven't already been covered

      if num_contibuting == 0:
        domain['low'] = None
        domain['high'] = None
        domain['covered'] = None
        continue

      domain['low'] = start_idx
      domain['high'] = start_idx + num_contibuting - 1
      domain['covered'] = num_covered

      fmax = self._seglen * width             # increment the frequencies that have been covered
      start_idx = start_idx + num_contibuting # increment the starting index for next time


  def _preferred_segment(self, f):
    # find a segment containing the desired frequency
    # prefer the highest reolution bins possible
    seg = None
    pow = len(self._segments)
    while seg is None:
      pow -= 1
      if pow < 0:
        break

      s = self._segments[pow]
      if s is None:
        continue

      _, width = s
      if f <= (width * self._seglen):
        seg = s
    
    if seg is None:
      raise ValueError('f outside domain of fft')
    return seg

  def track(self, pow):
    if pow >= self._seglen:
      raise ValueError('max pow is lg(n)')
    if self._segments[pow] is None:
      self._segments[pow] = FFTQuery.Segment([0] * self._seglen, self._base/(0x01 << pow), {'low': None, 'high': None, 'covered': None})
    # now recompute domains
    self._recompute_domains()
  
  def set(self, pow, data):
    if len(data) != self._seglen:
      raise ValueError('input data must have len n/2')
    seg = self._segments[pow]
    if seg is None:
      raise ValueError('no tracking for this pow')
    buf, _, _ = seg
    for idx, v in enumerate(data):
      buf[idx] = v


  def normalize(self):
    # takes absolute value
    # scales everything by a constant value such that max(self) is 1.0
    max = 0
    for idx, (_, s) in enumerate(self):
      a = abs(s)
      self._internal__setitem__(idx, a)
      if a > max:
        max = a
  
    # use max value to scale
    for idx, (_, s) in enumerate(self):
      self._internal__setitem__(idx, s/max)

  def interpolate(self, f):
    seg = self._preferred_segment(f)
    b, width, _ = seg
    idxl, idxh = FFTQuery._bookend_idxs(f, width, self._seglen)

    fl, fh = (idxl * width, idxh * width) # low and high bin freqs
    sl, sh = (b[idxl], b[idxh])           # low and high strengths
    delta = f - fl                        # distance from low freq

    return (f, delta * ((sh - sl)/(fh - fl)) + sl) # linear interp

  def closest(self, f):
    seg = self._preferred_segment(f)
    b, width = seg
    idx = FFTQuery._closest_idx(f, width, self._seglen)
    return (width * idx, b[idx])
  
  def floor(self, f):
    seg = self._preferred_segment(f)
    b, width = seg
    idx = FFTQuery._floor_idx(f, width)
    return (width * idx, b[idx])
  
  def strongest(self):
    a = [abs(s) for _, s in self]
    strongest = max(a)
    return self[a.index(strongest)]
    


def main():
  import numpy as np
  import time
  import math
  import matplotlib.pyplot as plt
  import asyncio
  import sys
  import pyaudio # http://people.csail.mit.edu/hubert/pyaudio/

  ############
  # user setup

  # samples = 64
  samples = 512

  sample_freq = 16000
  # sample_freq = 44100
  # sample_freq = 192000

  state = {
    'signal_freq': 6000
  }

  audio_source = 'generator' # 'microphone' or 'generator'

  # end user setup
  ################

  def signal(state, sample_freq):
    count = 0
    while True:
      yield math.sin(2*math.pi * state['signal_freq'] * (count / sample_freq))
      count += 1
      if count > sample_freq:
        count = 0

  def sample(s, n):
    for _ in range(n):
      # yield from s()
      yield next(s)
      
  def bin_resolution(sample_freq, samples):
    return sample_freq/samples

  pa = pyaudio.PyAudio()
  stream = pa.open(format=pyaudio.paInt16,
                  channels=1,
                  rate=sample_freq,
                  input=True,
                  frames_per_buffer=samples)

  dec = Decimator(samples)
  fftq = FFTQuery(samples, sample_freq)

  # configure decimator and fft query tracking to the same powers

  # dec.decimate(2)
  # fftq.track(2)  

  dec.decimate(7)
  fftq.track(7)

  s = signal(state, sample_freq)

  def got_stdin_data(q):
    data = sys.stdin.readline()
    freq = float(data)
    state['signal_freq'] = freq
    asyncio.create_task(q.put(freq))

  async def print_freq(q):
    while True:
      freq = await q.get()
      print(freq)

  async def analyze():
    while True:
      try:
        fig, (ax1, ax2) = plt.subplots(1, 2)

        if audio_source == 'generator':
          source = list(sample(s, samples))
        elif audio_source == 'microphone':
          audio = stream.read(samples, exception_on_overflow=False)
          source = []
          for i in range(len(audio)//2):
            source.append(int.from_bytes(audio[2*i:2*i+1], byteorder='little'))
        else:
          raise Exception('unknown audio source!')

        print(len(source))

        dec.feed(source)

        for pow, b in dec:
          ax2.plot(b)

          fs = sample_freq/(0x01 << pow)
          strengths = [abs(np.real(s)) for s in np.fft.fft(b)]
          freqs = np.fft.fftfreq(samples, 1/fs)
          z = filter(lambda z: z[0] >= 0, zip(freqs, strengths))
          freqs, strengths = zip(*z)

          fftq.set(pow, strengths)
        
        fftq.normalize()
        print(fftq.strongest())

        strengths = []
        freqs = []
        for strength, freq in fftq:
          strengths.append(strength)
          freqs.append(freq)
        
        ax1.plot(strengths, freqs, marker=None)
        plt.show(block=False)
        plt.pause(0.05)
        plt.close()

        await asyncio.sleep(samples/sample_freq)

        while True:
          try:
            pass
          except KeyboardInterrupt:
            break
          
      except KeyboardInterrupt:
        # quit()
        break

  print('decimator')
  print('')
  print('sample frequency: ', sample_freq/1000, ' kHz')
  print('')
  print('        decimator: ', end='')
  for idx, _ in enumerate(dec):
    print(idx, end=', ')
  print('')
  print('                   ---------')
  print('   division power: ', end='')
  for pow, b in dec:
    print(pow, end=', ')
  print('')
  print('   bin sizes (hz): ', end='')
  for pow, b in dec:
    print(bin_resolution(sample_freq/(0x01 << pow), len(b)), end=', ')
  print('')
  print('        max freqs: ', end='')
  for pow, b in dec:
    print((sample_freq/(0x01 << pow) /2), end=', ')
  print('\n')


  print(
  'the program will compute the strengths on each decimator\n'
  'and output the strongest detected frequency\n'
  )
  print('')
  print(
    'you may enter a new signal frequency on the commandline \n'
    'and press enter to see it take effect\n'
  )

  c = input('press enter to continue')

  loop = asyncio.get_event_loop()
  q = asyncio.Queue()
  loop.add_reader(sys.stdin, got_stdin_data, q)
  loop.create_task(print_freq(q))
  loop.create_task(analyze())

  try:
    loop.run_forever()
  except KeyboardInterrupt:
    quit()

if __name__ == '__main__':
  main()
