#!/usr/bin/env python3

from collections import deque

def highest_power_of_two(n):
  count = 0
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
    self._bufs = [None] * highest_power_of_two(n)
    self.decimate(0)

  def __getitem__(self, idx):
    return self._bufs[idx]
  
  def __iter__(self):
    return ((pow, b) for pow, b in enumerate(self._bufs) if b is not None)

  def decimate(self, pow):
    if (0x01 << pow) > (self._buflen):
      raise ValueError('decimator power division cannot exceed sample length')
    self._bufs[pow] = DecimationBuffer(self._buflen)

  def feed(self, data):
    for pow, buf in enumerate(self._bufs):
      if buf is not None:
        if pow == 0:
          buf.feed(data)
        else:
          skips = 0x01 << pow
          buf.feed(data[::skips])


def main():
  import numpy as np
  import time
  import math
  import matplotlib.pyplot as plt
  import asyncio
  import sys

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

  samples = 512
  sample_freq = 44000
  # sample_freq = 192000
  state = {
    'signal_freq': 6000
  }

  dec = Decimator(samples)
  # dec.decimate(2)
  # dec.decimate(4)
  dec.decimate(8)

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
      source = list(sample(s, samples))
      dec.feed(source)

      for pow, b in dec:
        fs = sample_freq/(0x01 << pow)
        
        fft = [abs(np.real(s)) for s in np.fft.fft(b)]
        freqs = np.fft.fftfreq(samples, 1/fs)
        z = filter(lambda z: z[0] > 0, zip(freqs, fft))
        freqs, fft = zip(*z)

        max_idx = fft.index(max(fft))
        print(f'{float(freqs[max_idx]):6.4f}', end=', ')
      
      print('')

      await asyncio.sleep(samples/sample_freq)

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
  'the program will compute the FFT on each decimator\n'
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
    pass

if __name__ == '__main__':
  main()
