import logging

# variable to control the logging level of the simulator
yarf_logging_level = logging.ERROR

formatter = logging.Formatter('[%(levelname)s]: %(name)s: %(message)s')

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)

""" Yarf Logger """
yarf_logger = logging.getLogger('yarf.Yarf')
yarf_logger.setLevel(yarf_logging_level)

yarf_logger.addHandler(stream_handler)

""" DC Analysis Logger """
dc_logger = logging.getLogger('Yarf.Analyses.DC')
dc_logger.setLevel(yarf_logging_level)

dc_logger.addHandler(stream_handler)

#file_handler = logging.FileHandler('DC.log')
#file_handler.setFormatter(formatter)
#dc_logger.addHandler(file_handler)

""" AC Analysis Logger """
ac_logger = logging.getLogger('Yarf.Analyses.AC')
ac_logger.setLevel(yarf_logging_level)

ac_logger.addHandler(stream_handler)

#file_handler = logging.FileHandler('AC.log')
#file_handler.setFormatter(formatter)
#ac_logger.addHandler(file_handler)

""" Transient Analysis Logger """
tr_logger = logging.getLogger('Yarf.Analyses.Transient')
tr_logger.setLevel(yarf_logging_level)

tr_logger.addHandler(stream_handler)

#file_handler = logging.FileHandler('AC.log')
#file_handler.setFormatter(formatter)
#ac_logger.addHandler(file_handler)
