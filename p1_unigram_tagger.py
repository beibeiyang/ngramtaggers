#! /usr/bin/python

__author__="Beibei Yang <byang1@cs.uml.edu>"

from collections import defaultdict

def replace_rare_words (oldtrainfile, newtrainfile, rarewords):
    """
    read from oldtrainfile, replace the rare words with _RARE_,
    and output to newtrainfile
    """
    fin = open(oldtrainfile, 'r')
    fout = open(newtrainfile, 'w')
    numline = 0
    numrare = 0
    for inline in fin:
        fields = inline.split(" ")
        if len(fields) == 2:
            word = fields[0].strip()
            if word in rarewords:
                fields[0] = '_RARE_'
                numrare += 1
                inline = ' '.join(fields)
        fout.write(inline)
        numline += 1
    fin.close()
    fout.close()
    print '%d lines written to %s\n' % (numline, newtrainfile)
    print 'Number of rare words: %d\n' % (numrare)
    
def get_rare_words (train_file):
    """
    Read from a file with lines like "1 WORDTAG O mind", 
    save the words if frequency count < 5
    """
    wordfreqs = defaultdict(int)
    rarewords = []
    f = open(train_file, 'r')
    for line in f:
        fields = line.split(" ")
        if len(fields) != 2:
            continue
        word = fields[0].strip()
        wordfreqs [word] += 1
    f.close()    
    #print wordfreqs
    for word, count in wordfreqs.iteritems():
        #print word, count
        if count < 5:
            # rare word
            rarewords.append(word)
    
    #print rarewords
    return rarewords


def rare_words_handler (train_file, new_train_file):
    """
    creates a new file new_train_file with all the rare words 
    in train_file replaced by _RARE_ 
    """
    rarewords = get_rare_words(train_file)
    replace_rare_words (train_file, new_train_file, rarewords)


def tag_counts (count_file):
    """
    Count(y) for each tag y 
    Output in the form of {'I-GENE': 41072, 'O': 345128})
    """
    tagcounts = defaultdict(int)
    f = open(count_file, 'r')
    for line in f:
        fields = line.split()
        if fields[1] != 'WORDTAG':
            continue
        count = int(fields[0])
        tag = fields[2]
        tagcounts[tag] += count  
    f.close()      
    return tagcounts
            
def word_tag_counts (count_file):
    """
    Count(y~>x) for each word tagged with y 
    Output in the form of {'aggression': [('O', 6)], 'hormone': [('I-GENE', 56), ('O', 37)],...}
    """
    wordtagcounts = defaultdict(list)
    f = open(count_file, 'r')
    for line in f:
        fields = line.split()
        if fields[1] != 'WORDTAG':
            continue
        count = int(fields[0])
        tag = fields[2]
        word = fields[3] 
        wordtagcounts[word].append((tag, count)) 
    f.close()      
    return wordtagcounts

def get_emission (count_file):
    """
    Calculates the emission parameters
    Output: {'localizes': ('O', 2.3179805753227787e-05), 'yellow': ('O', 2.0282330034074315e-05), ...}
    """
    emission_parameters = defaultdict(list)
    wordtagcounts = word_tag_counts (count_file)
    tagcounts = tag_counts (count_file)
    for word, tag_cnts in wordtagcounts.iteritems():
        emission_parameters [word] = (None, -1)
        for tag, cnt in tag_cnts:
            emission = float(cnt) / tagcounts[tag]
            if emission > emission_parameters [word][1]:
                emission_parameters[word] = (tag, emission)                
    return emission_parameters        
        
def pos_tagging (emission, corpus_file, output_file):
    fin = open (corpus_file, 'r')
    fout = open (output_file, 'w')
    
    for line in fin:
        word = line.strip()
        if len(word) < 1:
            fout.write (line)
            continue
        tag = ''
        if emission[word]:
            tag = emission[word][0]
        else:
            tag = emission['_RARE_'][0]  
        fout.write (' '.join([word, tag]) + '\n')
    fin.close()
    fout.close() 
        

if __name__ == "__main__":
    # prepare 'gene.train.RARE' for 'gene.counts.RARE'
    rare_words_handler('gene.train', 'gene.train.RARE')
   
    #emission = get_emission ('gene.counts.RARE')
    #pos_tagging (emission, 'gene.dev', 'gene_dev.p1.out')
    #pos_tagging (emission, 'gene.test', 'gene_test.p1.out')
#    