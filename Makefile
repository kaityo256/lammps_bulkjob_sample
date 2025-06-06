all:
	@echo "run make clean for cleaning"

.PHONY: clean

clean:
	rm -rf data/T*
	rm -f data/job1.sh
	rm -f data/job8.sh
