from sldb.common.mutations import threshold_mutations
import sldb.util.funcs as funcs


class Exporter(object):
    def __init__(self, session, rtype, rids, export_fields, selected_fields):
        self.session = session
        self.rtype = rtype
        self.rids = rids
        self.selected_fields = selected_fields
        self.export_fields = export_fields

    def _name_and_field(self, e):
        """Gets the name and field for export field entry ``e``.

        :param tuple e: Input field to split into name and field.

        :returns: Name and field for export field ``e``
        :rtype: (name, field)

        """
        if type(e) == str:
            return e, lambda r: getattr(r, e)
        return e

    def get_selected_data(self, seq, **overrides):
        """Gets the data specified by ``selected_fields`` for the sequence
        ``seq`` while overriding values in ``overrides`` if they exist.

        :param Sequence seq: The sequence from which to gather fields
        :param kwargs overrides: Fields to override

        :returns: A list of ``(name, value)`` tuples with the selected data
        :rtype: list

        """
        data = {}
        for field in self.export_fields:
            n, f = self._name_and_field(field)
            if n in self.selected_fields:
                try:
                    if n in overrides:
                        data[n] = overrides[n]
                    elif 'clone' not in n or (seq.clone is not None):
                        data[n] = f(seq)
                    else:
                        data[n] = ''
                except:
                    data[n] = ''
        return data

    def get_headers(self):
        headers = []
        for field in self.export_fields:
            n, f = self._name_and_field(field)
            if n in self.selected_fields:
                headers.append(n)
        return headers

    def get_base_query(self, model):
        return self.session.query(model).filter(
            getattr(
                model, '{}_id'.format(self.rtype)
            ).in_(self.rids),
        )
