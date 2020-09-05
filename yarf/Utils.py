import logging

# variables to control the logging levels
yarf_level = logging.ERROR
dc_level   = logging.ERROR
ac_level   = logging.ERROR
tr_level   = logging.ERROR

formatter = logging.Formatter('[%(levelname)s]: %(name)s: %(message)s')

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)

""" Yarf Logger """
yarf_logger = logging.getLogger('yarf.Yarf')
yarf_logger.setLevel(yarf_level)

yarf_logger.addHandler(stream_handler)

""" DC Analysis Logger """
dc_logger = logging.getLogger('Yarf.Analyses.DC')
dc_logger.setLevel(dc_level)

dc_logger.addHandler(stream_handler)

# file_handler = logging.FileHandler('DC.log')
# file_handler.setFormatter(formatter)
# dc_logger.addHandler(file_handler)

""" AC Analysis Logger """
ac_logger = logging.getLogger('Yarf.Analyses.AC')
ac_logger.setLevel(ac_level)

ac_logger.addHandler(stream_handler)

#file_handler = logging.FileHandler('AC.log')
#file_handler.setFormatter(formatter)
#ac_logger.addHandler(file_handler)

""" Transient Analysis Logger """
tr_logger = logging.getLogger('Yarf.Analyses.Transient')
tr_logger.setLevel(tr_level)

tr_logger.addHandler(stream_handler)

#file_handler = logging.FileHandler('AC.log')
#file_handler.setFormatter(formatter)
#ac_logger.addHandler(file_handler)
