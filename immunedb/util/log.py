import logging

logger = logging.getLogger('immunedb')
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(
    logging.Formatter('%(asctime)s [%(levelname)s] %(message)s',
                      datefmt='%Y-%m-%d %H:%M:%S')
)
logger.addHandler(handler)
