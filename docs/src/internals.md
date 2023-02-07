# Notes on optimisation

- We have three different assignment variables:
    - current
    - best
    - proposal
- We need best because we might accept a proposal which is worse than the curent value.
- This is dealt within `accept_reject_update!()`.