# ---------------------------------------------------------------------------
# EXPERIMENTS WITH NEW LARGER ALIGNMENT
# ---------------------------------------------------------------------------
ny_aligned_LARGE = list(SeqIO.parse('new_compiled_aligned.txt', "fasta")) 
ny_aligned_original = list(SeqIO.parse('new_york_compiled.txt', "fasta")) 
record = ny_aligned_LARGE[0]
new_records = []

for i, s in enumerate(str(record.seq)):
    if s == '-':
        print(i)

NEW_BAD_SEQIDS = BAD_SEQ_IDS.copy()
        
GAP_INDICES = [11074, 11075,  11076, 28253, 28254, 28255, 29897, 29898, 29899, 
                29900, 29901, 29902, 29903, 29904, 29905, 29906, 29907, 29908]
print(len(ny_aligned_LARGE))
for record in ny_aligned_LARGE[1:]:
    for index in GAP_INDICES:
        if record.seq[index] != '-':
            if record.id not in NEW_BAD_SEQIDS:
                NEW_BAD_SEQIDS.append(str(record.id))
                print(record.id.split("|")[-1])

for record in ny_aligned_LARGE[1:]:
    if record.id not in NEW_BAD_SEQIDS:
        record.id = record.description
        new_records.append(record)
print(len(NEW_BAD_SEQIDS))

location_dict = {}
time_dict = {}

for record in new_records[1:]:
    location = record.id.split('|')[-1]
    time = record.id.split('|')[-2].split('-')[1]
    if location in location_dict:
        location_dict[location]+=1
    else:
        location_dict[location] = 1
    if time in time_dict:
        time_dict[time]+=1
    else:
        time_dict[time] = 1

print(time_dict)
print(location_dict)

final_records = []
for record in ny_aligned_original:
    if record.id not in NEW_BAD_SEQIDS:
        final_records.append(record)

SeqIO.write(final_records, "new_compiled.txt", 'fasta')
