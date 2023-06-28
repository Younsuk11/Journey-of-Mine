# Improving Performance of PRIDE Algorithm through Advanced Methods

**reference**
> Tigunova, Anna, et al. "PRIDE: Predicting Relationships in Conversations." Proceedings of the 2021 Conference on Empirical Methods in Natural Language Processing. 2021.
---

Duration : 2023.4.4~2023.6.13

Team
* Department : Biomedical Engineering
* Younsuk Yeom, Chaehyeon Kim, Donghyeon Ki

## 0. Task : Relationship Extracting Tasks
the task of identifying and classifying relationships between entities in text. 

## 1. Idea
**Limits in PRIDE**
> "Due to the 512 input lenght limitation of BERT"
* splits input sequence of utterances into chunks and runs BERT for each chunk
* it is difficult to say it as "conversational context"
* do not make good use of BERT in identifying relationship between chunks
* assumption that each chunk has maximal possible length

**Solution of ours**
1. Randomly suffle the order of utterances
   * exchange information between words
2. Create overlapped chunks
   * takes full advantage of BERt & understand "smooth" context
  
## 2. Contribution
1. Experiments for Two methods
2. Comparing each method's performances for each relationship label
3. Observing which method helps to classify specific relationship
4. Examining the effectiveness of PRIDE's hierarchical representation of conversation
   
