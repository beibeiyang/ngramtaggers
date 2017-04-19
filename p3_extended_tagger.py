#! /usr/bin/python

__author__ = "Beibei Yang <byang1@cs.uml.edu>"

from collections import defaultdict
import re
import sys

################################################################################
## Part 3 of assignment 1
## Step 1: 
##     Execute rare_words_handler to generate gene.train.RARE.4
## Step 2: 
##     python count_freqs.py gene.train.RARE.4 > gene.counts.RARE.4
## Step 3: 
##     Execute pos_tag to generate gene_dev.p3.out or gene_test.p3.out
## Step 4: To evaluate, execute
##     python eval_gene_tagger.py gene.key gene_dev.p3.out
################################################################################

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
    for word, count in wordfreqs.iteritems():
        #print word, count
        if count < 5:
            # rare word
            rarewords.append(word)
    return rarewords

def get_rare_class (rareword):
    if re.search('[0-9]+', rareword):
        return '_NUMERIC_'
    if re.search('^[A-Z]+$', rareword):
        return '_ALLCAPS_'
    if re.search('.*[A-Z]$', rareword):
        return '_LASTCAP_'
    return '_RARE_'

def set_rare_classes (rarewords):
    rarewords_c = defaultdict(lambda: '')
    for word in rarewords:
        cl = ''
        if word in rarewords:
            cl = get_rare_class (word)

        rarewords_c [word] = cl            
    return rarewords_c

def replace_rare_words (oldtrainfile, newtrainfile, rarewords_classes):
    """
    read from oldtrainfile, replace the rare words with _RARE_ or other rare tags,
    and output to newtrainfile
    """
    fin = open(oldtrainfile, 'r')
    fout = open(newtrainfile, 'w')
    numline = 0
    rareclass_counter = defaultdict(int)
    for inline in fin:
        fields = inline.split(" ")
        if len(fields) == 2:
            word = fields[0].strip()
            
            rareclass = rarewords_classes [word]
            if rareclass:    # word is in rarewords_classes
                # replace the word with the rareclass
                fields[0] = rareclass
                inline = ' '.join(fields)
                rareclass_counter[rareclass] += 1
        fout.write(inline)
        numline += 1
    fin.close()
    fout.close()
    print '%d lines written to %s\n' % (numline, newtrainfile)
    print str(rareclass_counter), "\n"
    
    
def rare_words_handler (train_file, new_train_file):
    """
    creates a new file new_train_file with all the rare words 
    in train_file replaced by _RARE_ 
    """
    rarewords = get_rare_words(train_file)
    rarewords = set_rare_classes (rarewords)
    print rarewords
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
        fields = line.split(" ")
        if fields[1] != 'WORDTAG':
            continue
        count = int(fields[0].strip())
        tag = fields[2].strip()
        word = fields[3].strip()
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
        #emission_parameters [word] = defaultdict(float)
        emission_parameters [word] = defaultdict(float)
        for tag, cnt in tag_cnts:
            emission = float(cnt) / tagcounts[tag]
            emission_parameters[word][tag] = emission               
    return emission_parameters 

def get_ngram_count (count_file, n=3):
    ngramcounts = defaultdict(list)
    assert n in xrange(1,4), "n = 1,2,3"
    f = open(count_file, 'r')
    ngram_tag = '-'.join([str(n), 'GRAM'])
    for line in f:
        fields = line.split()
        if fields[1] != ngram_tag:
            continue
        count = int(fields[0])
        ngram = ' '.join(fields[2:]).strip()
        ngramcounts[ngram] = count
    f.close()      
    return ngramcounts

def get_trigram_conditional (y1, y2, y3, trigram_counts, bigram_counts):
    """
    Computes q(y_{i}|y_{i-2}, y_{i-1}).
    This corresponds to the first bullet of part 2
    """
    return float(trigram_counts[ ' '.join([y1, y2, y3]) ]) / bigram_counts[ ' '.join([y1, y2]) ]


def viterbi ( X, count_file, trigram_counts, bigram_counts, emission, debug=False):
    n = len(X)
    X = [None] + X # counter starts at 1
    S = {}
    S[-1] = ['*']
    S[0] = ['*']
    for i in xrange(n):
        S[i+1] = ['O', 'I-GENE']
    
    pi = {}
    bp = {}
    pi[(0, '*', '*')] = 1.00
    emission = get_emission (count_file)
    if debug:
        print emission
    
    for k in xrange(1, n+1): # 1..n
        for u in S[k-1]:
            for v in S[k]:
                #pi = get_pi_max (k, u, v, S, X, emission)
                if debug:
                    print "k = %d \t U = %s \t V = %s" % (k, u, v)
                _pi = []
                for w in S[k-2]:
                    if debug:
                        print "For W = %s Calculating Pi[%d, %s, %s] * q(%s|%s,%s) * e(%s | %s)" % (w, k-1,w,u, v,w,u, X[k-1],v)

                    p = pi[(k-1,w,u)]   
                    q = get_trigram_conditional(w,u,v,trigram_counts,bigram_counts)
                    
                    e = 0
                    if emission[X[k]]: # list is not empty
                        e = emission[X[k]][v]
                    else:
                        # X_{k} is an unseen word, i.e. rare word
                        rareclass = get_rare_class (X[k])
                        e = emission[rareclass][v]
                        if debug:
                            print '  RARE word: %s \t Class: %s' % (X[k], rareclass)
                            
                    if debug:
                        print "  Pi[%d, %s, %s] = %f" % (k-1,w,u, p)
                        print "  q(%s|%s, %s) = %f" % (v,w,u, q)
                        print "  e(%s | %s) = %f" % (X[k],v, e)
                    _pi.append ( p * q * e )
                pi[(k,u,v)] = max(_pi)
                if debug:
                    print "  Probabilities: %s" % (str(_pi))
                    print "  pi[(%d, %s, %s)] = %s" % (k,u,v, pi[(k,u,v)]) 
                    print "  S[k-2]:", str(S[k-2]), '\t _pi.index( pi[(k,u,v)] ):', str(_pi.index( pi[(k,u,v)]))  
                                                                                     
                bp[(k,u,v)] = S[k-2][_pi.index( pi[(k,u,v)] )]
                if debug:
                    print "  bp[(%d, %s, %s)] = %s" % (k,u,v, bp[(k,u,v)]) 
                
#    print len(S), len(X)
#    print pi, bp

    max_pi, max_u, max_v = -1, None, None            
#    for u in S[1]:
#        for v in S[2]:
#            _pi = pi[(n, u, v)] * get_trigram_conditional(u,v,'STOP',trigram_counts,bigram_counts)
#            if _pi > max_pi:
#                max_pi = _pi
#                max_u, max_v = u, v
                
    for (_n, u, v) in pi:
        if _n is n:
            _pi = pi[(n, u, v)] * get_trigram_conditional(u,v,'STOP',trigram_counts,bigram_counts)
            if _pi >= max_pi:
                max_pi = _pi
                max_u, max_v = u, v

    y = [max_u, max_v]
    for k in xrange(n-2, 0, -1): # (n-2)..1
        _y = y[:]
        y = [bp[(k+2, _y[0], _y[1])]] + _y 
        #y = [bp[(k+2, y[0], y[1])]] + y

    return y
 
    
 
def simple_corpus_iterator(corpus_file):
    """
    Get an iterator object over the corpus file. The elements of the
    iterator contain (word, ne_tag) tuples. Blank lines, indicating
    sentence boundaries return (None, None).
    """
    l = corpus_file.readline()
    while l:
        line = l.strip()
        if line: # Nonempty line
            yield line
        else: # Empty line
            yield None                        
        l = corpus_file.readline()


def sentence_iterator(corpus_iterator):
    """
    Return an iterator object that yields one sentence at a time.
    Sentences are represented as lists of (word, ne_tag) tuples.
    """
    current_sentence = [] #Buffer for the current sentence
    for l in corpus_iterator:        
            if l is None:
                if current_sentence:  #Reached the end of a sentence
                    yield current_sentence
                    current_sentence = [] #Reset buffer
                else: # Got empty input stream
                    sys.stderr.write("WARNING: Got empty input file/stream.\n")
                    raise StopIteration
            else:
                current_sentence.append(l) #Add token to the buffer

    if current_sentence: # If the last line was blank, we're done
        yield current_sentence  #Otherwise when there is no more token
                                # in the stream return the last sentence.
  


def pos_tag (corpus_file, count_file, output_file):
    trigram_counts = get_ngram_count (count_file, 3)
    bigram_counts = get_ngram_count (count_file, 2)
    emission = get_emission(count_file)
    
    sent_iterator = sentence_iterator(simple_corpus_iterator (file (corpus_file, 'r')) )
    
    fout = open (output_file, 'w')

    for sent in sent_iterator:
        Y = viterbi ( sent, count_file, trigram_counts, bigram_counts, emission)
        
        for i in range(len(sent)):
            fout.write (' '.join([sent[i], Y[i]]) + '\n')
        
        fout.write ('\n')
    
    fout.close()
 
     
if __name__ == "__main__":
    # prepare 'gene.train.RARE' for 'gene.counts.RARE'
    #rare_words_handler('gene.train', 'gene.train.RARE.4')
    
    # part 3 of the assignment
    #pos_tag ( 'gene.test.s', 'gene.counts.RARE.4', 'gene_dev.p3.out')
    #pos_tag ( 'gene.dev', 'gene.counts.RARE.4', 'gene_dev.p3.out')
    pos_tag ( 'gene.test', 'gene.counts.RARE.4', 'gene_test.p3.out')
