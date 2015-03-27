import datetime
import json

from sldb.common.models import ModificationLog


class LoggingException(Exception):
    pass


def make_clone_mod(new_clone, old_v_name=None, old_v_seq=None, gaps=None):
    if old_v_name is None and gaps is None:
        raise LoggingException('No modifications specified')
    if old_v_name is not None and old_v_seq is None:
        raise LoggingException('Must specify old sequence with old_v_name')

    info = {
        'clone_id': new_clone.id,
    }
    if old_v_name != new_clone.v_gene:
        info['old_v_name'] = old_v_name
        info['old_v_seq'] = old_v_seq
        info['new_v_name'] = new_clone.v_gene
        info['new_v_seq'] = new_clone.group.germline

    if gaps is not None:
        info['gaps'] = gaps

    return make_mod('clone_modify', info=info)


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
