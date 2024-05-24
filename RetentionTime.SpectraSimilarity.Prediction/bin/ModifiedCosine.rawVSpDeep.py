from optparse import OptionParser
from matchms.importing import load_from_mgf
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms import calculate_scores
from matchms.similarity import CosineGreedy
from matchms.similarity import ModifiedCosine

def main():
    # Read spectrums from a MGF formatted file, for other formats see https://matchms.readthedocs.io/en/latest/api/matchms.importing.html
    # Apply default filter to standardize ion mode, correct charge and more.
    # Default filter is fully explained at https://matchms.readthedocs.io/en/latest/api/matchms.filtering.html .
    # Apply filters to clean and enhance each spectrum[raw]
    file = load_from_mgf(options.raw)
    spectrums = []
    for spectrum in file:
        spectrum = default_filters(spectrum)
        # Scale peak intensities to maximum of 1
        spectrum = normalize_intensities(spectrum)
        spectrums.append(spectrum)

    # Apply filters to clean and enhance each spectrum[pDeep2 no normalize_intensities]
    query_spectrums_file=load_from_mgf(options.predict)
    query_spectrums = []
    for query_spectrum in query_spectrums_file:
        query_spectrum = default_filters(query_spectrum)
        # Scale peak intensities to maximum of 1
        # query_spectrum = normalize_intensities(query_spectrum)
        query_spectrums.append(query_spectrum)

    # Calculate Cosine similarity scores between all spectrums
    # For other similarity score methods see https://matchms.readthedocs.io/en/latest/api/matchms.similarity.html .
    # Because references and queries are here the same spectra, we can set is_symmetric=True
    # scores = calculate_scores(references=spectrums,
    #                           queries=query_spectrums,
    #                           similarity_function=CosineGreedy())

    similarity_measure = ModifiedCosine(tolerance=0.05)

    # scores = calculate_scores(spectrums, query_spectrums, similarity_measure)

    result = open(options.outfile, 'w')
    print('{}\t{}\t{}'.format( "raw_spectral", "pDeep2_pred","similarity_score"), file = result)
    for i in range(len(spectrums)):
        for j in range(len(query_spectrums)):
            if (spectrums[i].metadata["title"].split('#')[0] == query_spectrums[j].metadata["title"].split('#')[0]) and (spectrums[i].metadata["title"].split('#')[1] == query_spectrums[j].metadata["title"].split('#')[1]) and (spectrums[i].metadata["title"].split('#')[2] == query_spectrums[j].metadata["title"].split('#')[2]):
                scores = calculate_scores([spectrums[i]], [query_spectrums[j]], similarity_measure)
                print('{}\t{}\t{}'.format(spectrums[i].metadata["title"], query_spectrums[j].metadata["title"], float((scores.scores[0][0]['score']))), file = result)
                break

if __name__ == "__main__":
    parser = OptionParser(prog='coolpython',
                            usage="%prog [OPTION]... FILE...",
                            description="**************",)
    parser.add_option("-i", "--raw",dest="raw", help="")
    parser.add_option("-p", "--predict",dest="predict", help="")
    parser.add_option("-o", "--outfile", dest="outfile", help="")
    (options, args) = parser.parse_args()
    main()