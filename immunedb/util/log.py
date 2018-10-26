import coloredlogs
import logging

logger = logging.getLogger('immunedb')
colors = coloredlogs.DEFAULT_FIELD_STYLES
colors['levelname']['color'] = 'white'
coloredlogs.install(level='DEBUG', logger=logger,
                    fmt='%(asctime)s [%(levelname)s] %(message)s',
                    field_styles=colors)
