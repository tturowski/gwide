.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/gwide/master/docs/_static/images/title.png
   :alt: [logo]
   :align: center
   :width: 400px

gwide
==========

(the documentation is under construction, come back later for more)


.. code:: python

    from gwide import (SeparatingDigestionsProblem, list_common_enzymes,
                           LADDERS, load_genbank)


    # DEFINE SEQUENCES AND ENZYME SET (6-CUTTERS WITH >3 COMMERCIAL PROVIDERS)
    enzymes = list_common_enzymes(site_length=(6,), min_suppliers=3)
    sequences = [
        load_genbank(genbank_file_path, name=f)
        for genbank_file_path in some_llist_of_files)
    ]

    # SELECT THE BEST DIGESTION PAIRS (AT MOST 1 ENZYME PER DIGESTION)
    problem = SeparatingDigestionsProblem(enzymes=enzymes,
                                          ladder=LADDERS['100_to_4k'],
                                          sequences=sequences,
                                          max_enzymes_per_digestion=1)
    score, selected_digestions = problem.select_digestions(max_digestions=2)

    # GENERATE A FIGURE OF THE BAND PATTERNS

    problem.plot_digestions(
        selected_digestions,
        patterns_props={'label_fontdict': {'rotation': 35}},
        target_file="separating_digestions.png"
    )

    problem.plot_distances_map(digestions=selected_digestions,
                               target_file="separating_digestions_distances.png")

Usage: Construct validation or identification from experimental data
---


Contribute
----------

gwide is an open-source library originally written at the
Edinburgh Genome Foundry by Zulko_. It is released on Github_ under the MIT
licence (¢ Edinburgh Genome Foundry), everyone is welcome to contribute.

.. _Zulko: https://github.com/Zulko/
.. _Github: https://github.com/EdinburghGenomeFoundry/gwide
.. _PYPI: https://pypi.python.org/pypi/gwide


.. raw:: html

       <a href="https://twitter.com/share" class="twitter-share-button"
       data-text="gwide - Digestion enzyme selection with Python" data-size="large" data-hashtags="Bioprinting">Tweet
       </a>
       <script>!function(d,s,id){var js,fjs=d.getElementsByTagName(s)[0],p=/^http:/.test(d.location)?'http':'https';
       if(!d.getElementById(id)){js=d.createElement(s);js.id=id;js.src=p+'://platform.twitter.com/widgets.js';
       fjs.parentNode.insertBefore(js,fjs);}}(document, 'script', 'twitter-wjs');
       </script>
       <iframe src="http://ghbtns.com/github-btn.html?user=Edinburgh-Genome-Foundry&repo=gwide&type=watch&count=true&size=large"
       allowtransparency="true" frameborder="0" scrolling="0" width="152px" height="30px" margin-bottom="30px"></iframe>


.. toctree::
    :hidden:
    :maxdepth: 3

    self

.. toctree::
    :hidden:
    :caption: Reference
    :maxdepth: 3

    ref


.. _Zulko: https://github.com/Zulko/
.. _Github: https://github.com/EdinburghGenomeFoundry/gwide
.. _PYPI: https://pypi.python.org/pypi/gwide
