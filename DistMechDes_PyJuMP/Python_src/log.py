"""
Customized logging for writing to log files and consoles
Handlet .setLevel determines the logging level to write to file or console
Logging levels are as follows:

Level   Numerical Value
CRITICAL    50
ERROR   40
WARNING 30
INFO    20
DEBUG   10
NOTSET  10
"""
import logging
import os

log = logging.getLogger('Simulate Horizontal Investment Coordination for Transmission Expansion Planning in Decentralized mode')
log.setLevel(logging.DEBUG)

logfile = os.path.join(os.getcwd(), "Simulate Horizontal Investment Coordination for Transmission Expansion Planning in Decentralized mode.log")

file_handler = logging.FileHandler(filename=logfile, mode='a')
file_formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
file_handler.setFormatter(file_formatter)
file_handler.setLevel(logging.INFO)

console_handler = logging.StreamHandler()
console_formatter = logging.Formatter('%(name)-12s %(levelname)-8s %(message)s')
console_handler.setFormatter(console_formatter)
console_handler.setLevel(logging.INFO)

log.addHandler(file_handler)
log.addHandler(console_handler)
