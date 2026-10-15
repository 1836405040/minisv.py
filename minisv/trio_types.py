"""Type-aware candidate normalization and explicit pending structural review."""


def normalize_record(fields):
    fields = list(fields)
    info = dict(x.split('=', 1) for x in fields[7].split(';') if '=' in x)
    typ = info.get('SVTYPE', '').upper()
    info['SVTYPE'] = typ
    if typ == 'TRA':
        # Preserve the original representation; a TRA alias alone does not
        # establish a valid paired BND junction or orientation.
        return fields, info, None, 'UNSUPPORTED_TRA_REPRESENTATION'
    if typ not in ('INS', 'DEL', 'INV', 'DUP', 'BND'):
        return fields, info, None, 'UNSUPPORTED_SVTYPE'
    if typ == 'BND':
        return fields, info, None, None
    try:
        length = int(info.get('SVLEN', 0))
        if not length and typ in ('DEL', 'DUP', 'INV'):
            length = int(info.get('END', fields[1])) - int(fields[1])
            if length <= 0:
                return fields, info, None, 'UNRESOLVED_SVLEN'
        if not length:
            return fields, info, None, 'UNRESOLVED_SVLEN'
    except ValueError:
        return fields, info, None, 'INVALID_SVLEN_OR_END'
    if typ == 'DEL':
        length = -abs(length)
    info['SVLEN'] = str(length)
    flags = [x for x in fields[7].split(';') if x and '=' not in x]
    fields[7] = ';'.join(flags + [f'{k}={v}' for k, v in info.items()])
    return fields, info, abs(length), None


def partition_for_review(source_path, destination, pending_path):
    """Retain unsupported records with provenance; return supported record count."""
    import gzip
    import json
    opener = gzip.open if str(source_path).endswith('.gz') else open
    count = 0
    with opener(source_path, 'rt') as source, open(destination, 'w') as target, open(pending_path, 'a') as pending:
        for line in source:
            if line.startswith('#'):
                target.write(line)
                continue
            if not line.strip():
                continue
            fields, info, length, error = normalize_record(line.rstrip().split('\t'))
            if error or info['SVTYPE'] not in ('INS', 'DEL'):
                pending.write(json.dumps(dict(id=fields[2], type=info['SVTYPE'],
                    source=str(source_path), status=error or 'PENDING_TYPE_SPECIFIC_REVIEW',
                    original_record=line.rstrip(),
                    reason='No validated type-specific final review model')) + '\n')
            else:
                target.write('\t'.join(fields) + '\n')
                count += 1
    return count
