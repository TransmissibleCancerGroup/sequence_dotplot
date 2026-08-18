import numpy as np 
import matplotlib.pyplot as plt

number = {'A':0, 'T':1, 'G':2, 'C':3}

def read_to_numb(F):
    with open(f"{F}.txt") as f:
        text = f.read().replace("\n", "")
    number_list = [number[i] for i in list(text)]
    number_mat = np.array(number_list)
    return number_mat

#get read sequences from text files
seq1 = read_to_numb('sequence1')
seq2 = read_to_numb('sequence2')

n1 = len(seq1)
n2 = len(seq2)

print('sequence 1 length:',n1)
print('sequence 2 length:',n2)

#load options
options_dict = {}
labels_dict = {}
with open('options.txt') as options: 
    for line in options:
        line = line.rstrip('\n')
        line = line.rstrip(' ')
        if len(line.split(':')) > 1:
            label, text = line.split(':')
            if text and all(c.isdigit() or c == ',' for c in text):
                coords = [int(x) for x in text.split(",")]
                options_dict[label] = coords
            else:
                labels_dict[label] = text
options.close()

if 'window_size' in options_dict:
    w = int(options_dict['window_size'][0]) #window size
    print('windows size:',w)
else: 
    assert 'Please input a window size'

X = []
Y = []
X_rc = []
Y_rc = []

# ------ main sequence matchup ------
#
#my original code, which chatgpt has optimised. Keeping it for reference as it is easier to understand 
#
# for i in range(n1-w):
#     if i%100 == 0:
#         print(i)
#     for j in range(i,n2-w):
#         if np.array_equal(seq1[i:(i+w)],seq2[j:(j+w)]):
#             X.append(i)
#             Y.append(j)

# plt.scatter(X,Y, s=1, c='blue')
# plt.scatter(Y,X, s=1, c='blue')
# plt.show()


# ---- optimised: precompute window hashes (credits ChatGPT) ----
hash1 = [hash(seq1[i:i+w].tobytes()) for i in range(n1 - w + 1)]
hash2 = [hash(seq2[i:i+w].tobytes()) for i in range(n2 - w + 1)]

# map hash -> positions in seq2
hash_to_pos = {}
for j, h in enumerate(hash2):
    hash_to_pos.setdefault(h, []).append(j)

# ---- main loop (structure preserved) ----

print(f'Progress out of {n1} bases:')
print('Forward:')
for i in range(n1 - w):
    if i % 10000 == 0:
        print(i)
    h = hash1[i]
    if h in hash_to_pos:
        for j in hash_to_pos[h]:
            # optional collision safety check
            if np.array_equal(seq1[i:i+w], seq2[j:j+w]):
                X.append(i)
                Y.append(j)

# complement mapping:
# A(0)<->T(1), G(2)<->C(3)
complement = np.array([1, 0, 3, 2])

# reverse complement of seq2
seq2_rc = complement[seq2[::-1]]

# hashes of reverse complement
hash2_rc = [hash(seq2_rc[i:i+w].tobytes()) for i in range(n2 - w + 1)]

# map hash -> positions in reverse complement
hash_to_pos_rc = {}
for j, h in enumerate(hash2_rc):
    hash_to_pos_rc.setdefault(h, []).append(j)

print('Reverse:')
# ---- reverse complement matches ----
for i in range(n1 - w):
    if i % 10000 == 0:
        print(i)

    h = hash1[i]

    if h in hash_to_pos_rc:
        for j in hash_to_pos_rc[h]:

            # collision safety check
            if np.array_equal(seq1[i:i+w], seq2_rc[j:j+w]):

                X_rc.append(i)

                # convert RC coordinate back to original seq2 coordinate
                Y_rc.append(n2 - w - j)

print('Plotting...')

#label gene regions: 

if 'sequence_1_start' in options_dict:
    start_pos_1 = options_dict['sequence_1_start'][0] 
else:
    start_pos_1 = 0
    print('No options for sequence_1_start found, starting at 0')

if 'sequence_2_start' in options_dict:
    start_pos_2 = options_dict['sequence_2_start'][0] 
else:
    start_pos_2 = 0
    print('No options for sequence_2_start found, starting at 0')

#plotting 
X = [x + start_pos_1 for x in X]
Y = [x + start_pos_2 for x in Y]
X_rc = [x + start_pos_1 for x in X_rc]
Y_rc = [x + start_pos_2 for x in Y_rc]

# forward matches
plt.scatter(X, Y, s=0.1, c='blue', label='forward')

# reverse complement matches
plt.scatter(X_rc, Y_rc, s=0.1, c='red', label='reverse complement')

#label regions on plot:
cmap = plt.get_cmap('hsv') 
# ^ change this for a different colour scheme (google 'matplotlip colourmaps')
# depending on the chosen colourmap, may also need to change the range of the list below
colours = [cmap(x) for x in np.linspace(0,1,len(options_dict)-3, endpoint=False)]
c = 0
for label, region in options_dict.items():
    if not label in ['window_size','sequence_1_start','sequence_2_start']:
        if region[0] > start_pos_1 and start_pos_1 + n1 > region[1]: #condition for gene to be on the scale
            plt.axvspan(region[0], region[1], color=colours[c], alpha=0.3, label=label)
            c += 1

plt.ticklabel_format(style='plain', axis='both', useOffset=False)

plt.title(labels_dict['title'])
plt.xlabel(labels_dict['x_label'])
plt.ylabel(labels_dict['y_label'])

plt.legend()
plt.show()