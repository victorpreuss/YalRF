import logging

# variables to control the logging levels
yalrf_level = logging.ERROR
dc_level   = logging.ERROR
ac_level   = logging.ERROR
tr_level   = logging.ERROR

formatter = logging.Formatter('[%(levelname)s]: %(name)s: %(message)s')

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)

""" YalRF Logger """
yalrf_logger = logging.getLogger('yalrf.YalRF')
yalrf_logger.setLevel(yalrf_level)

yalrf_logger.addHandler(stream_handler)

""" DC Analysis Logger """
dc_logger = logging.getLogger('YalRF.Analyses.DC')
dc_logger.setLevel(dc_level)

dc_logger.addHandler(stream_handler)

# file_handler = logging.FileHandler('DC.log')
# file_handler.setFormatter(formatter)
# dc_logger.addHandler(file_handler)

""" AC Analysis Logger """
ac_logger = logging.getLogger('YalRF.Analyses.AC')
ac_logger.setLevel(ac_level)

ac_logger.addHandler(stream_handler)

#file_handler = logging.FileHandler('AC.log')
#file_handler.setFormatter(formatter)
#ac_logger.addHandler(file_handler)

""" Transient Analysis Logger """
tr_logger = logging.getLogger('YalRF.Analyses.Transient')
tr_logger.setLevel(tr_level)

tr_logger.addHandler(stream_handler)

#file_handler = logging.FileHandler('AC.log')
#file_handler.setFormatter(formatter)
#ac_logger.addHandler(file_handler)
