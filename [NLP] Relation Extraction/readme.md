# Improving Performance of PRIDE Algorithm through Advanced Methods

**reference**
> Tigunova, Anna, et al. "PRIDE: Predicting Relationships in Conversations." Proceedings of the 2021 Conference on Empirical Methods in Natural Language Processing. 2021.
---

Team
* Department : Biomedical Engineering
* Younsuk Yeom, Chaehyeon Kim, Donghyeon Ki

## 0. Task : Relationship Extracting Tasks
the task of identifying and classifying relationships between entities in text. 

## 1. Idea
**Limits in PRIDE**
"Due to the 512 input lenght limitation of BERT"
* splits input sequence of utterances into chunks and runs BERT for each chunk
* it is difficult to say it as "conversational context"
* do not make good use of BERT in identifying relationship between chunks
* assumption that each chunk has maximal possible length
