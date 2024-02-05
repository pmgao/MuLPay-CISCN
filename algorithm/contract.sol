pragma solidity ^0.8.0;

contract MultiSigWallet {
    uint public numConfirmationsRequired;
    address[] public owners;
    mapping(address => bool) public isOwner;
    uint public transactionCount;
    mapping(uint => Transaction) public transactions;
    mapping(uint => mapping(address => bool)) public isConfirmed;

    struct Transaction {
        address to;
        uint value;
        bytes data;
        bool executed;
        uint numConfirmations;
    }

    event Deposit(address indexed sender, uint value);
    event Submission(uint indexed transactionId);
    event Confirmation(address indexed sender, uint indexed transactionId);
    event Execution(uint indexed transactionId);
    event ExecutionFailure(uint indexed transactionId);

    constructor(address[] memory _owners, uint _numConfirmationsRequired) {
        require(_owners.length > 0, "Owners required");
        require(
           _numConfirmationsRequired > 0 && _numConfirmationsRequired <= _owners.length,
            "Invalid number of required confirmations"
        );

        for (uint i = 0; i < _owners.length; i++) {
            address owner = _owners[i];
            require(owner != address(0), "Invalid owner");
            require(!isOwner[owner], "Duplicate owner");

            isOwner[owner] = true;
            owners.push(owner);
        }

        numConfirmationsRequired = _numConfirmationsRequired;
    }

    receive() external payable {
        emit Deposit(msg.sender, msg.value);
    }

    modifier onlyOwner() {
        require(isOwner[msg.sender], "Not an owner");
        _;
    }

    modifier transactionExists(uint _transactionId) {
        require(transactions[_transactionId].to != address(0), "Transaction does not exist");
        _;
    }

    modifier confirmed(uint _transactionId) {
        require(isConfirmed[_transactionId][msg.sender], "Transaction not confirmed");
        _;
    }

    function submitTransaction(address _to, uint _value, bytes memory _data) public onlyOwner {
        uint transactionId = transactionCount++;

        transactions[transactionId] = Transaction({
            to: _to,
            value: _value,
            data: _data,
            executed: false,
            numConfirmations: 0
        });

        emit Submission(transactionId);
        confirmTransaction(transactionId);
    }

   function confirmTransaction(uint _transactionId) public onlyOwner transactionExists(_transactionId) {
        Transaction storage transaction = transactions[_transactionId];

        require(!isConfirmed[_transactionId][msg.sender], "Transaction already confirmed");

        isConfirmed[_transactionId][msg.sender] = true;
        transaction.numConfirmations++;

        emit Confirmation(msg.sender, _transactionId);

        if (transaction.numConfirmations >= numConfirmationsRequired) {
            executeTransaction(_transactionId);
        }
    }

    function executeTransaction(uint _transactionId) public onlyOwner transactionExists(_transactionId) confirmed(_transactionId) {
        Transaction storage transaction = transactions[_transactionId];

        require(!transaction.executed, "Transaction already executed");

        transaction.executed = true;

        (bool success, ) = transaction.to.call{value: transaction.value}(transaction.data);
        if (success) {
            emit Execution(_transactionId);
        } else {
            emit ExecutionFailure(_transactionId);
            transaction.executed = false;
        }
    }
}