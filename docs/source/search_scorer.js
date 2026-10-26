/*
 * Search scorer for the PolyzyMD docs (`html_search_scorer` in conf.py).
 *
 * Sphinx scores a full-text hit by the best single query word: a page with
 * any one word in any heading scores the same as a page that has every word
 * in its title. This scorer keeps the Sphinx weights and adds to them:
 *
 *   +5 for each query word in the page title,
 *   +3 for each query word in a section heading of the page,
 *   +3 for pages in get_started/, tutorials/ and how_to/, which answer tasks,
 *   -5 for api/ and contributor_guide/ pages and the generated indexes.
 *
 * Index entries (the `{index}` directive) still score 50 to 100 and stay on
 * top. scripts/check_docs_search.py measures the effect on real queries.
 */
var Scorer = {
  // Sphinx defaults (sphinx/themes/basic/static/searchtools.js).
  objNameMatch: 11,
  objPartialMatch: 6,
  objPrio: { 0: 15, 1: 5, 2: -5 },
  objPrioDefault: 0,
  title: 15,
  partialTitle: 7,
  term: 5,
  partialTerm: 2,

  score: (result) => {
    const [docname, , , , score] = result;
    let bonus = 0;
    if (/^(get_started|tutorials|how_to)\//.test(docname)) bonus += 3;
    if (/^(api|contributor_guide)\/|^(genindex|py-modindex)$/.test(docname)) bonus -= 5;

    const index = typeof Search !== "undefined" ? Search._index : null;
    const query = new URLSearchParams(window.location.search).get("q");
    if (!index || !query) return score + bonus;
    const file = index.docnames.indexOf(docname);
    if (file < 0) return score + bonus;

    const stemmer = new Stemmer();
    const stem = (word) => stemmer.stemWord(word.toLowerCase());
    const words = (text) =>
      text.toLowerCase().split(/[^\p{L}\p{N}]+/u).filter((w) => w.length > 2 && !stopwords.includes(w));
    const titleWords = new Set(words(index.titles[file]).map(stem));
    const inHeading = (w) => {
      const files = index.titleterms[w];
      return files !== undefined && (files === file || (files.includes && files.includes(file)));
    };
    for (const w of new Set(words(query).map(stem))) {
      if (titleWords.has(w)) bonus += 5;
      if (inHeading(w)) bonus += 3;
    }
    return score + bonus;
  },
};
