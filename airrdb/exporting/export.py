class Exporter(object):
    def __init__(self, session, rtype, rids, export_fields, selected_fields):
        self.session = session
        self.rtype = rtype
        self.rids = rids
        self.selected_fields = selected_fields
        self.export_fields = export_fields

    def get_selected_data(self, seq):
        """Gets the data specified by ``selected_fields`` for the sequence
        ``seq`` while overriding values in ``overrides`` if they exist.

        :param Sequence seq: The sequence from which to gather fields

        :returns: A list of ``(name, value)`` tuples with the selected data
        :rtype: list

        """
        return {
            name: accessor(seq) for name, accessor in
            self.export_fields.iteritems() if name in self.selected_fields
        }
