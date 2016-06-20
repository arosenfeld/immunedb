import datetime
import json

from immunedb.common.models import ModificationLog


class LoggingException(Exception):
    pass


def make_mod(action_type, info, session=None, commit=False):
    ml = ModificationLog(
        datetime=datetime.datetime.utcnow(),
        action_type=action_type,
        info=json.dumps(info))
    if session is None:
        if commit:
            raise LoggingException('Specified `commit` without `session`')
        return ml

    session.add(ml)
    if commit:
        session.commit()
    return ml
