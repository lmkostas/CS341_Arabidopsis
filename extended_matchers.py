import os
import re
import warnings
# Travis will not import the PorterStemmer
if 'CI' not in os.environ:
    try:
        from nltk.stem.porter import PorterStemmer
    except ImportError:
        warnings.warn("nltk not installed- some default functionality may be absent.")
    '''
    try:
        from nltk.stem.wordnet import WordNetLemmatizer
    except ImportError:
        warnings.warn("nltk not installed- some default functionality may be absent.")
    '''


class Matcher(object):
    """
    Applies a function f : c -> {True,False} to a generator of candidates,
    returning only candidates _c_ s.t. _f(c) == True_,
    where f can be compositionally defined.
    """
    def __init__(self, *children, **opts):
        self.children           = children
        self.opts               = opts
        self.longest_match_only = self.opts.get('longest_match_only', True)
        self.init()
        self._check_opts()

    def init(self):
        pass

    def _check_opts(self):
        """
        Checks for unsupported opts, throws error if found
        NOTE: Must be called _after_ init()
        """
        for opt in self.opts.keys():
            if not self.__dict__.has_key(opt):
                raise Exception("Unsupported option: %s" % opt)

    def _f(self, c):
        """The internal (non-composed) version of filter function f"""
        return True

    def f(self, c):
        """
        The recursively composed version of filter function f
        By default, returns logical **conjunction** of operator and single child operator
        """
        if len(self.children) == 0:
            return self._f(c)
        elif len(self.children) == 1:
            return self._f(c) and self.children[0].f(c)
        else:
            raise Exception("%s does not support more than one child Matcher" % self.__name__)

    def _is_subspan(self, c, span):
        """Tests if candidate c is subspan of span, where span is defined specific to candidate type"""
        return False

    def _get_span(self, c):
        """Gets a tuple that identifies a span for the specific candidate class that c belongs to"""
        return c

    def apply(self, candidates):
        """
        Apply the Matcher to a **generator** of candidates
        Optionally only takes the longest match (NOTE: assumes this is the *first* match)
        """
        seen_spans = set()
        for c in candidates:
            if self.f(c) and (not self.longest_match_only or not any([self._is_subspan(c, s) for s in seen_spans])):
                if self.longest_match_only:
                    seen_spans.add(self._get_span(c))
                yield c


WORDS = 'words'

class NgramMatcher(Matcher):
    """Matcher base class for Ngram objects"""
    def _is_subspan(self, c, span):
        """Tests if candidate c is subspan of span, where span is defined specific to candidate type"""
        return c.char_start >= span[0] and c.char_end <= span[1]

    def _get_span(self, c):
        """Gets a tuple that identifies a span for the specific candidate class that c belongs to"""
        return (c.char_start, c.char_end)


class DictionaryMatch(NgramMatcher):
    """Selects candidate Ngrams that match against a given list d"""
    def init(self):
        self.ignore_case = self.opts.get('ignore_case', True)
        self.attrib      = self.opts.get('attrib', WORDS)
        self.reverse     = self.opts.get('reverse', False)
        try:
            self.d = frozenset(w.lower() if self.ignore_case else w for w in self.opts['d'])
        except KeyError:
            raise Exception("Please supply a dictionary (list of phrases) d as d=d.")

        # Optionally use a stemmer, preprocess the dictionary
        # Note that user can provide *an object having a stem() method*
        self.blacklist = self.opts.get('blacklist', None)
        self.stemmer = self.opts.get('stemmer', None)
        if self.stemmer is not None:
            if self.stemmer == 'porter':
                self.stemmer = PorterStemmer()
            if self.blacklist is not None:
                self.d = frozenset(self._stem(w) for w in list(self.d) if self._stem(w) not in self.blacklist and len(w)>2)
            else:
                self.d = frozenset(self._stem(w) for w in list(self.d))

    def _stem(self, w):
        """Apply stemmer, handling encoding errors"""
        try:
            if len(w.split()) == 1:
                return self.stemmer.stem(w)
            else:
                new_w = []
                for word in w.split():
                    new_w.append(self._stem(word))
                return ' '.join(new_w)
        except:# UnicodeDecodeError:
            return w

    def _f(self, c):
        p = c.get_attrib_span(self.attrib)
        p = p.lower() if self.ignore_case else p
        p = self._stem(p) if self.stemmer is not None else p        
        return (not self.reverse) if p in self.d else self.reverse

class LambdaFunctionMatch(NgramMatcher):
    """Selects candidate Ngrams that match against a given list d"""
    def init(self):
        self.ignore_case = self.opts.get('ignore_case', True)
        self.attrib      = self.opts.get('attrib', WORDS)
        try:
            self.func = self.opts['func']
        except KeyError:
            raise Exception("Please supply a dictionary (list of phrases) d as d=d.")
    
    def _f(self, c):
        """The internal (non-composed) version of filter function f"""
        return self.func(c)


class Union(NgramMatcher):
    """Takes the union of candidate sets returned by child operators"""
    def f(self, c):
       for child in self.children:
           if child.f(c) > 0:
               return True
       return False


class Concat(NgramMatcher):
    """
    Selects candidates which are the concatenation of adjacent matches from child operators
    NOTE: Currently slices on **word index** and considers concatenation along these divisions only
    """
    def init(self):
        self.permutations   = self.opts.get('permutations', False)
        self.left_required  = self.opts.get('left_required', True)
        self.right_required = self.opts.get('right_required', True)
        self.ignore_sep     = self.opts.get('ignore_sep', True)
        self.sep            = self.opts.get('sep', " ")

    def f(self, c):
        if len(self.children) != 2:
            raise ValueError("Concat takes two child Matcher objects as arguments.")
        if not self.left_required and self.children[1].f(c):
            return True
        if not self.right_required and self.children[0].f(c):
            return True

        # Iterate over candidate splits **at the word boundaries**
        for wsplit in range(c.get_word_start()+1, c.get_word_end()+1):
            csplit = c.word_to_char_index(wsplit) - c.char_start  # NOTE the switch to **candidate-relative** char index

            # Optionally check for specific separator
            if self.ignore_sep or c.get_span()[csplit-1] == self.sep:
                c1 = c[:csplit-len(self.sep)]
                c2 = c[csplit:]
                if self.children[0].f(c1) and self.children[1].f(c2):
                    return True
                if self.permutations and self.children[1].f(c1) and self.children[0].f(c2):
                    return True
        return False


class SlotFillMatch(NgramMatcher):
    """Matches a slot fill pattern of matchers _at the character level_"""
    def init(self):
        self.attrib = self.opts.get('attrib', WORDS)
        try:
            self.pattern = self.opts['pattern']
        except KeyError:
            raise Exception("Please supply a slot-fill pattern p as pattern=p.")

        # Parse slot fill pattern
        split        = re.split(r'\{(\d+)\}', self.pattern)
        self._ops    = map(int, split[1::2])
        self._splits = split[::2]

        # NOTE: Must have non-null splits!!
        if any([len(s) == 0 for s in self._splits[1:-1]]):
            raise ValueError("SlotFillMatch must have non-empty split patterns to function currently.")

        # Check for correct number of child matchers / slots
        if len(self.children) != len(set(self._ops)):
            raise ValueError("Number of provided matchers (%s) != number of slots (%s)." \
                    % (len(self.children), len(set(self._ops))))

    def f(self, c):

        # First, filter candidates by matching splits pattern
        m = re.match(r'(.+)'.join(self._splits) + r'$', c.get_attrib_span(self.attrib))
        #m = re.match(r'(.+)'.join(self._splits), c.get_attrib_span(self.attrib))
        if m is None:
            return False

        # Then, recursively apply matchers
        
        for i,op in enumerate(self._ops):
            if self.children[op].f(c[m.start(i+1):m.end(i+1)]) == 0:
                return False
        return True


class RegexMatch(NgramMatcher):
    """Base regex class- does not specify specific semantics of *what* is being matched yet"""
    def init(self):
        try:
            self.rgx = self.opts['rgx']
        except KeyError:
            raise Exception("Please supply a regular expression string r as rgx=r.")
        self.ignore_case = self.opts.get('ignore_case', True)
        self.attrib      = self.opts.get('attrib', WORDS)
        self.sep         = self.opts.get('sep', " ")
        self.blacklist  = self.opts.get('blacklist', None)

        # Compile regex matcher
        # NOTE: Enforce full span matching by ensuring that regex ends with $!
        self.rgx = self.rgx if self.rgx.endswith('$') else self.rgx + r'$'
        self.r = re.compile(self.rgx, flags=re.I if self.ignore_case else 0)

    def _f(self, c):
        raise NotImplementedError()


class RegexMatchSpan(RegexMatch):
    """Matches regex pattern on **full concatenated span**"""
    def _f(self, c):
        return True if self.r.match(c.get_attrib_span(self.attrib, sep=self.sep)) is not None else False


class RegexMatchEach(RegexMatch):
    """Matches regex pattern on **each token**"""
    def _f(self, c):
        tokens = c.get_attrib_tokens(self.attrib)
        return True if tokens and all([self.r.match(t) is not None for t in tokens]) else False    


class PersonMatcher(RegexMatchEach):
    """
    Matches Spans that are the names of people, as identified by CoreNLP.

    A convenience class for setting up a RegexMatchEach to match spans
    for which each token was tagged as a person.
    """
    def __init__(self, *children, **kwargs):
        kwargs['attrib'] = 'ner_tags'
        kwargs['rgx'] = 'PERSON'
        super(PersonMatcher, self).__init__(*children, **kwargs)


class LocationMatcher(RegexMatchEach):
    """
    Matches Spans that are the names of locations, as identified by CoreNLP.

    A convenience class for setting up a RegexMatchEach to match spans
    for which each token was tagged as a location.
    """
    def __init__(self, *children, **kwargs):
        kwargs['attrib'] = 'ner_tags'
        kwargs['rgx'] = 'LOCATION'
        super(LocationMatcher, self).__init__(*children, **kwargs)


class OrganizationMatcher(RegexMatchEach):
    """
    Matches Spans that are the names of organizations, as identified by CoreNLP.

    A convenience class for setting up a RegexMatchEach to match spans
    for which each token was tagged as an organization.
    """
    def __init__(self, *children, **kwargs):
        kwargs['attrib'] = 'ner_tags'
        kwargs['rgx'] = 'ORGANIZATION'
        super(OrganizationMatcher, self).__init__(*children, **kwargs)


class DateMatcher(RegexMatchEach):
    """
    Matches Spans that are dates, as identified by CoreNLP.

    A convenience class for setting up a RegexMatchEach to match spans
    for which each token was tagged as a date.
    """
    def __init__(self, *children, **kwargs):
        kwargs['attrib'] = 'ner_tags'
        kwargs['rgx'] = 'DATE'
        super(DateMatcher, self).__init__(*children, **kwargs)


class NumberMatcher(RegexMatchEach):
    """
    Matches Spans that are numbers, as identified by CoreNLP.

    A convenience class for setting up a RegexMatchEach to match spans
    for which each token was tagged as a number.
    """
    def __init__(self, *children, **kwargs):
        kwargs['attrib'] = 'ner_tags'
        kwargs['rgx'] = 'NUMBER'
        super(NumberMatcher, self).__init__(*children, **kwargs)


class MiscMatcher(RegexMatchEach):
    """
    Matches Spans that are miscellaneous named entities, as identified by CoreNLP.

    A convenience class for setting up a RegexMatchEach to match spans
    for which each token was tagged as miscellaneous.
    """
    def __init__(self, *children, **kwargs):
        kwargs['attrib'] = 'ner_tags'
        kwargs['rgx'] = 'MISC'
        super(MiscMatcher, self).__init__(*children, **kwargs)

class Sequence(NgramMatcher):
    """Matches a slot fill pattern of matchers _at the character level_"""
    def init(self):
        self.attrib = self.opts.get('attrib', WORDS)
        self.sep            = self.opts.get('sep', ' ')
        self.links            = self.opts.get('links', None)
        self.required = self.opts.get('required', [1 for _ in range(len(self.children))])
        self.links  =self.opts.get('links', None)

        # Check for correct number of child matchers / slots
        if len(self.children) < 1:
            raise ValueError("Number of provided matchers is 0.")

    def f(self, c):
        # First, filter candidates by matching splits pattern
        match_instance = [(r+1)%2 for r in self.required]
        m = re.match(r'([\w-]+)', c.get_attrib_span(self.attrib))
        #print m, c.get_attrib_span(self.attrib)
        if m is None:
            return False
        # Then, recursively apply matchers
        rgx = r'\S+'
        if self.links is not None:
            rgx = r'([A-Za-z0-9-]+)'
        for match in re.finditer(rgx, c.get_attrib_span(self.attrib)):
            match_found = False
            #print match, c.get_attrib_span(self.attrib)
            split = c.get_attrib_span(self.attrib)[match.start(0):match.end(0)]
            for i in range(len(self.children)):
                if self.children[i].f(c[match.start(0):match.end(0)]) != 0:
                    match_found = True
                    match_instance[i] = 1
                    break
            if not match_found:
                #return False
                if self.links:
                    #if (split != 'and' and split != 'or' and split != 'and/or'): return False
                    for l in range(len(self.links)):
                        if self.links[l].f(c[match.start(0):match.end(0)]) != 0:
                            match_found = True
                            break
                    if not match_found:
                        return False
                else:
                    return False
        if sum(match_instance) != len(match_instance): return False
        return True

class Contains(NgramMatcher):
    """Matches a slot fill pattern of matchers _at the character level_"""
    def init(self):
        try:
            self.rgxs = self.opts['rgxs']
        except KeyError:
            raise Exception("Please supply a regular expression string r as rgx=r.")
        self.ignore_case = self.opts.get('ignore_case', True)
        self.attrib      = self.opts.get('attrib', [WORDS]*len(self.rgxs))
        self.sep         = self.opts.get('sep', " ")
        self.type =self.opts.get('type', 'and')

        # Compile regex matcher
        # NOTE: Enforce full span matching by ensuring that regex ends with $!
        self.rgxs = [rgx if rgx.endswith('$') else rgx + r'$' for rgx in self.rgxs]
        self.r = [re.compile(rgx, flags=re.I if self.ignore_case else 0) for rgx in self.rgxs]

        # Check for correct number of child matchers / slots
        if len(self.children) != 1:
            raise ValueError("Number of provided matchers is not 1.")

    def f(self, c):
        # First, filter candidates by matching splits pattern
        for child in self.children:
            if child.f(c) <= 0: return False

        if self.type =='or':
            for i, rgx in enumerate(self.r):
                if re.search(rgx, c.get_attrib_span(self.attrib[i], sep=self.sep)) is not None:
                    return True

            return False
        else:
            for i, rgx in enumerate(self.r):
                if re.search(rgx, c.get_attrib_span(self.attrib[i], sep=self.sep)) is None:
                    return False 

            return True

"""
class Contains(NgramMatcher):
    '''Matches a slot fill pattern of matchers _at the character level_'''
    def init(self):
        '''
        try:
            self.rgxs = self.opts['rgxs']
        except KeyError:
            raise Exception("Please supply a regular expression string r as rgx=r.")
        '''
        try:
            self.query = self.opts['query']
        except KeyError:
            raise Exception("Please supply a query q as query=q.")
        try:
            self.type = self.opts['type']
        except KeyError:
            raise Exception("Please supply a type for query t as type=t.")
        self.ignore_case = self.opts.get('ignore_case', True)
        self.attrib      = self.opts.get('attrib', WORDS)
        self.sep         = self.opts.get('sep', " ")
        self.reverse    =self.opts.get('reverse', False)

        # Compile regex matcher
        # NOTE: Enforce full span matching by ensuring that regex ends with $!
        if self.type == 'rgx':
            self.rgxs = [rgx if rgx.endswith('$') else rgx + r'$' for rgx in self.query]
            self.r = [re.compile(rgx, flags=re.I if self.ignore_case else 0) for rgx in self.rgxs]

        # Check for correct number of child matchers / slots
        if len(self.children) > 1:
            raise ValueError("Number of provided matchers is greater than 1.")

    def f(self, c):
        # First, filter candidates by matching splits pattern
        if self.children:
            for child in self.children:
                if child.f(c) <= 0: return False
        if not self.reverse:
            if self.type == 'rgx':
                for rgx in self.r:
                    if re.search(rgx, c.get_attrib_span(self.attrib, sep=self.sep)) is None:
                        return False 
                return True
            elif self.type == 'words':
                for d in self.query:
                    for match in re.finditer(r'\w+', c.get_attrib_span('lemmas')):
                        if re.match(d, c[match.start(1):match.end(1)], flags = re.IGNORECASE) is None:
                            return False
            return True
        else:
            if self.type == 'rgx':
                for rgx in self.r:
                    if re.search(rgx, c.get_attrib_span(self.attrib, sep=self.sep)) is not None:
                        return False 
            elif self.type =='words':
                for match in re.finditer(r'\w+', c.get_attrib_span('lemmas')):
                    for d in self.query:
                        if re.match(d, c.get_attrib_span(self.attrib)[match.start(0):match.end(0)], flags = re.IGNORECASE) is not None:
                            return False
            return True
"""
class Intersection(NgramMatcher):
    """Takes the intersection of candidate sets returned by child operators"""
    def f(self, c):
       for child in self.children:
           if child.f(c) <= 0:
               return False
       return True

class MergeMatcher(NgramMatcher):
    """
    Looks at overlapping or adjacent subsets of the span.
    Runs two child matchers on overlapping subsets, and returns true if there is some overlapping
    where each of the children return true for their part.
    Compare to: Concat only returns true if they are next to each other, but doesn't check overlaps.
    """
    def init(self):
        self.permutations   = self.opts.get('permutations', False)

    def f(self, c):
        if len(self.children) != 2:
            raise ValueError("MergeMatcher takes two child Matcher objects as arguments.")
        # Iterate over candidate splits **at the word boundaries**
        # Split 1
        for wsplit in range(c.get_word_start()+1, c.get_word_end()+1):
            # Split 2
            for wsplit2 in range(c.get_word_start(), wsplit+1):

                csplit = c.word_to_char_index(wsplit) - c.char_start  # NOTE the switch to **candidate-relative** char index
                csplit2 = c.word_to_char_index(wsplit2) - c.char_start

                c1 = c[:csplit]
                c2 = c[csplit2:]
                if self.children[0].f(c1) and self.children[1].f(c2):
                    return True
                if self.permutations and self.children[1].f(c1) and self.children[0].f(c2):
                    return True
        return False

class Seq(NgramMatcher):
    def init(self):
        self.ignore_case = self.opts.get('ignore_case', True)
        self.attrib      = self.opts.get('attrib', WORDS)
        self.reverse     = self.opts.get('reverse', False)
        try:
            self.d = list(set(w.lower() if self.ignore_case else w for w in self.opts['d']))
            self.d.sort(key=len, reverse=True)
        except KeyError:
            raise Exception("Please supply a dictionary (list of phrases) d as d=d.")

        # Optionally use a stemmer, preprocess the dictionary
        # Note that user can provide *an object having a stem() method*
        self.stemmer = self.opts.get('stemmer', None)
        if self.stemmer is not None:
            if self.stemmer == 'porter':
                self.stemmer = PorterStemmer()
            self.d = frozenset(self._stem(w) for w in list(self.d))

    def _stem(self, w):
        """Apply stemmer, handling encoding errors"""
        try:
            return self.stemmer.stem(w)
        except UnicodeDecodeError:
            return w

    def _f(self, c):
        p = c.get_attrib_span(self.attrib)
        p = p.lower() if self.ignore_case else p
        p = self._stem(p) if self.stemmer is not None else p
        while len(p) > 0: 
            for w in self.d:
                if len(p) > len(w):
                    if p[:len(w)] == w:
                        p = p[len(w):]
                elif len(p) == len(w) and p == w:
                    p = ''
        return len(p) == 0
